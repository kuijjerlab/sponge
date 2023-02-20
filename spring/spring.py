import pandas as pd
import numpy as np
import time

from pyjaspar import jaspardb

from uniprot_api import *

from math import log2


def plogp(x):
    """
    Returns x*log2(x) for a number, handles the 0 case properly.
    """
    if x == 0:
        return 0
    else:
        return x*log2(x)


def calculate_ic(motif):
    """
    Calculates the information content for a given motif, assumes equal ACGT distribution.
    """
    df = pd.DataFrame(motif.pwm)
    df['IC'] = df.apply(lambda x: 2 + sum([plogp(x[y]) for y in ['A', 'C', 'G', 'T']]), axis=1)
    
    return df['IC'].sum()


def adjust_name(i):
    return i[:-2] + i[-2:].lower()


class Spring:
    def __init__(self):
        pass


    def select_tfs(self):
        # Database object
        jdb_obj = jaspardb()
        print ('Using:', jdb_obj.release)

        # All vertebrate motifs
        motifs = jdb_obj.fetch_motifs(collection='CORE', tax_group='vertebrates', all_versions=True)
        print ('All motif versions:', len(motifs))
        print ('Motif base IDs:', len(set([i.base_id for i in motifs])))

        # Select latest, preferring human ones
        latest = {}
        for i in motifs:
            if i.base_id not in latest:
                latest[i.base_id] = [i.matrix_id, i.species]
            else:
                # Replace with newer version if the new one is human or the old one isn't       
                if ('9606' in i.species) or ('9606' not in latest[i.base_id][1]):
                    # This could be added to the logical condition above but this is more readable
                    if int(i.matrix_id[-1]) > int(latest[i.base_id][0][-1]):
                        latest[i.base_id] = [i.matrix_id, i.species]
        motifs_latest = [i for i in motifs if i.matrix_id == latest[i.base_id][0]]
        
        # Keep only one motif per TF
        # Consider dropping this requirement maybe
        tf_to_motif = {}
        for i in motifs_latest:
            if i.name not in tf_to_motif:
                tf_to_motif[i.name] = {i.matrix_id: calculate_ic(i)}
            else:
                tf_to_motif[i.name][i.matrix_id] = calculate_ic(i)
        self.tf_to_motif = tf_to_motif
        motifs_unique = [i for i in motifs_latest if tf_to_motif[i.name][i.matrix_id] == max(tf_to_motif[i.name].values())]
        print ('Unique motifs:', len(motifs_unique))

        # Drop heterodimers
        motifs_nohd = [i for i in motifs_unique if ':' not in i.name]
        print ('Motifs without heterodimers:', len(motifs_nohd))

        self.motifs = motifs_nohd


    def find_human_homologs(self, path_to_homologene='homologene.tsv'):
        # Get the non-human motifs
        non_human_motifs = [i for i in self.motifs if '9606' not in i.species]
        print ('Non-human motifs:', len(non_human_motifs))

        # Read the homologene database
        hg_df = pd.read_csv(path_to_homologene, sep='\t', header=None, 
                    names=['HG Group ID', 'TaxID', 'Gene ID', 
                           'Gene Symbol', 'Protein GI', 'Protein Accession'])

        # Get the non-human motif names
        non_human_motif_names = [i.name for i in non_human_motifs]
        # Compare against homologene
        found_names = hg_df[hg_df['Gene Symbol'].isin([adjust_name(i) for i in non_human_motif_names])]['Gene Symbol'].unique()
        # Find the missing ones
        missing = set([adjust_name(i) for i in non_human_motif_names]) - set(found_names)
        print ('Names missing from the homologene database:')
        for i in [(i.name,i.acc) for i in non_human_motifs if i.name in missing]:
            print (i[0], i[1])

        # Get the missing IDs from Uniprot API
        print ('Retrieving results from UniProt...')
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="RefSeq_Protein", ids=[i.acc[0] for i in non_human_motifs if i.name in missing]
        )
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)

        # Convert the results into a DataFrame
        mapping = pd.DataFrame(results['results'])
        mapping.columns = ['Uniprot', 'Accession']

        # Create a DataFrame for matching missing entries
        hg_df[hg_df['Protein Accession'].isin(mapping['Accession'])]
        missing_df = pd.DataFrame([(i.name, i.acc[0]) for i in non_human_motifs if i.name in missing],
            columns = ['Gene Symbol', 'Uniprot'])        
        matching_df = missing_df.join(mapping.set_index('Uniprot'), on='Uniprot').join(
            hg_df.set_index('Protein Accession'), on='Accession', rsuffix='_HG')

        def corresponding_id(name):
            values = hg_df[hg_df['Gene Symbol'] == name]['HG Group ID'].values
            if len(values) == 0:
                values = matching_df[matching_df['Gene Symbol'] == name]['HG Group ID'].values
            return values

        # Create a DataFrame of corresponding names
        corr_df = pd.DataFrame(non_human_motif_names, columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(adjust_name)
        corr_df['Group ID'] = corr_df['Adjusted Name'].apply(corresponding_id)
        corr_df['Group ID'] = corr_df['Group ID'].apply(lambda x: x[0] if len(x) > 0 else np.nan)
        corr_df['Human Name'] = corr_df['Group ID'].apply(lambda x: hg_df[(hg_df['HG Group ID'] == x) & 
                                                        (hg_df['TaxID'] == 9606)]['Gene Symbol'].values)
        corr_df['Human Name'] = corr_df['Human Name'].apply(lambda x: x[0] if len(x) > 0 else '')
        corr_df['Trivial'] = corr_df['Original Name'].apply(lambda x: x.upper()) == corr_df['Human Name']

        # Find duplicates
        duplicated = corr_df[corr_df['Human Name'].duplicated(keep=False) & (corr_df['Human Name'] != '')].copy()
        print ('Duplicate names:')
        print (duplicated)

        # Calculate the information content for duplicates
        duplicated['IC'] = duplicated['Original Name'].apply(lambda x: max(self.tf_to_motif[x].values()))
        # Keep the highest IC amongst the duplicates
        to_drop = duplicated['Original Name'][duplicated.sort_values('IC').duplicated('Human Name', keep='last')]

        # Exlude the IDs which are already present among the human ones
        human_motif_names = [i.name for i in self.motifs if '9606' in i.species]
        corr_df['Duplicate'] = corr_df['Human Name'].isin(human_motif_names)

        # Perform the final filtering - discard all duplicates and TFs without homologs
        corr_df_final = corr_df[(corr_df['Duplicate'] == False) & (corr_df['Human Name'] != '') & 
                        (corr_df['Original Name'].isin(to_drop) == False)]

        # The mapping of original to human names and the matrix IDs to be kept
        animal_to_human = {i: j for i,j in zip(corr_df_final['Original Name'], corr_df_final['Human Name'])}
        print ('Final number of IDs which will be replaced by human homologs:', len(animal_to_human))
        vert_list_nohd = [i.matrix_id for i in self.motifs if (i.name in human_motif_names
                          or i.name in animal_to_human)]
        print ('Final number of total matrix IDs:', len(vert_list_nohd))

        self.matrix_ids = vert_list_nohd
        self.animal_to_human = animal_to_human


    