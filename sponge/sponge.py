### Imports ###
import pandas as pd
import numpy as np
import bioframe
import time
import os
import requests

from pyjaspar import jaspardb

from multiprocessing import Pool

from math import *

from sponge.uniprot_api import *
from sponge.functions import *

### Class definition ###
class Sponge:
    def __init__(
        self,
        temp_folder = '.sponge_temp/'
    ):

        self.temp_folder = temp_folder


    def retrieve_file(
        self,
        description,
        temp_folder
    ):
        
        # TODO: Promoter file needs modifications
        # TODO: No reason this should be an instance function
        
        if description not in file_df.index:
            print (f'File description not recognised: {description}')
            return
        if not os.path.exists(temp_folder):
            os.mkdir(temp_folder)
        file_name = file_df.loc[description, 'name']
        file_path = os.path.join(temp_folder, file_name)
        if os.path.exists(file_path):
            print (f'Using cached file {file_path}')
            print ()
        else:
            print (f'File {file_name} not found in directory {temp_folder}')
            key = None
            positive = ['y', 'yes', 'hell yeah']
            negative = ['n', 'no', 'nope']
            while key is None or key.lower() not in positive + negative:
                if key is not None:
                    print (f'Input not recognised: {key}')
                print ('Do you want to download it? Y/N')
                key = input()
                print (key)
            print ()
            if key.lower() in negative:
                return None
            r = requests.get(file_df.loc[description, 'url'])
            mode = file_df.loc[description, 'mode']
            if mode == 'w':
                data = r.text
            else:
                data = r.content
            print (f'Downloading data into {file_path}...')
            open(file_path, mode).write(data)
            print ()     

        return file_path


    def select_tfs(
        self
    ):

        # Database object
        jdb_obj = jaspardb()
        self.release = jdb_obj.release
        print ('Using:', self.release)

        # All vertebrate motifs
        motifs = jdb_obj.fetch_motifs(collection='CORE', 
            tax_group='vertebrates', all_versions=True)
        print ()
        print ('All motif versions:', len(motifs))
        print ('Motif base IDs:', len(set([i.base_id for i in motifs])))

        # Select latest, preferring human ones
        latest = {}
        for i in motifs:
            if i.base_id not in latest:
                latest[i.base_id] = [i.matrix_id, i.species]
            else:
                # Replace with newer version if the new one is human or the old 
                # one isn't       
                if (('9606' in i.species) or 
                    ('9606' not in latest[i.base_id][1])):
                    # This could be added to the logical condition above but 
                    # this is more readable
                    if int(i.matrix_id[-1]) > int(latest[i.base_id][0][-1]):
                        latest[i.base_id] = [i.matrix_id, i.species]
        motifs_latest = [i for i in motifs if 
            i.matrix_id == latest[i.base_id][0]]
        
        # Keep only one motif per TF
        # Consider dropping this requirement maybe
        tf_to_motif = {}
        for i in motifs_latest:
            if i.name not in tf_to_motif:
                tf_to_motif[i.name] = {i.matrix_id: calculate_ic(i)}
            else:
                tf_to_motif[i.name][i.matrix_id] = calculate_ic(i)
        self.tf_to_motif = tf_to_motif
        motifs_unique = [i for i in motifs_latest if 
            (tf_to_motif[i.name][i.matrix_id] == 
            max(tf_to_motif[i.name].values()))]
        print ('Unique motifs:', len(motifs_unique))

        # Drop heterodimers
        motifs_nohd = [i for i in motifs_unique if ':' not in i.name]
        print ('Motifs without heterodimers:', len(motifs_nohd))

        self.motifs = motifs_nohd


    def find_human_homologs(
        self, 
        homologene_file=None
    ):
        
        if homologene_file is None:
            homologene_file = self.retrieve_file('homologene', self.temp_folder)
            if homologene_file is None:
                print ('Unable to find or retrieve the homologene file, exiting')
                return

        # Get the non-human motifs
        non_human_motifs = [i for i in self.motifs if '9606' not in i.species]
        print ('Non-human motifs:', len(non_human_motifs))

        # Read the homologene database
        hg_df = pd.read_csv(homologene_file, sep='\t', header=None, 
                    names=['HG Group ID', 'TaxID', 'Gene ID', 
                           'Gene Symbol', 'Protein GI', 'Protein Accession'])

        # Get the non-human motif names
        non_human_motif_names = [i.name for i in non_human_motifs]
        # Compare against homologene
        found_names = hg_df[hg_df['Gene Symbol'].isin([adjust_name(i) for i in 
            non_human_motif_names])]['Gene Symbol'].unique()
        # Find the missing ones
        missing = (set([adjust_name(i) for i in non_human_motif_names]) - 
            set(found_names))
        print ()
        print ('Names missing from the homologene database:')
        for i in [(i.name,i.acc) for i in non_human_motifs if 
            i.name in missing]:
            print (i[0], *i[1])

        # Get the missing IDs from Uniprot API
        print ()
        print ('Retrieving matches from UniProt...')
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="RefSeq_Protein", ids=[i.acc[0] 
            for i in non_human_motifs if i.name in missing]
        )
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)

        # Convert the results into a DataFrame
        mapping = pd.DataFrame(results['results'])
        mapping.columns = ['Uniprot', 'Accession']

        # Create a DataFrame for matching missing entries
        hg_df[hg_df['Protein Accession'].isin(mapping['Accession'])]
        missing_df = pd.DataFrame([(i.name, i.acc[0]) for i in non_human_motifs 
            if i.name in missing], columns = ['Gene Symbol', 'Uniprot'])        
        matching_df = missing_df.join(mapping.set_index('Uniprot'), 
            on='Uniprot').join(hg_df.set_index('Protein Accession'), 
            on='Accession', rsuffix='_HG')

        def corresponding_id(name):
            values = hg_df[hg_df['Gene Symbol'] == name]['HG Group ID'].values
            if len(values) == 0:
                matches = matching_df[matching_df['Gene Symbol'] == name]
                values = matches['HG Group ID'].values
            return values

        # Create a DataFrame of corresponding names
        corr_df = pd.DataFrame(non_human_motif_names, 
            columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(adjust_name)
        corr_df['Group ID'] = corr_df['Adjusted Name'].apply(corresponding_id)
        corr_df['Group ID'] = corr_df['Group ID'].apply(lambda x: 
            x[0] if len(x) > 0 else np.nan)
        corr_df['Human Name'] = corr_df['Group ID'].apply(lambda x: 
            hg_df[(hg_df['HG Group ID'] == x) & 
            (hg_df['TaxID'] == 9606)]['Gene Symbol'].values)
        corr_df['Human Name'] = corr_df['Human Name'].apply(lambda x: 
            x[0] if len(x) > 0 else '')
        corr_df['Trivial'] = corr_df['Original Name'].apply(lambda x: 
            x.upper()) == corr_df['Human Name']

        # Find duplicates
        duplicated = corr_df[corr_df['Human Name'].duplicated(keep=False) & 
            (corr_df['Human Name'] != '')].copy()
        to_print = duplicated.groupby('Human Name')['Original Name'].unique(
            ).apply(lambda x: ' '.join(x))
        print ()
        print ('Duplicate names:')
        for i in to_print.index:
            print (f'{i}:', to_print.loc[i])

        # Calculate the information content for duplicates
        duplicated['IC'] = duplicated['Original Name'].apply(lambda x: 
            max(self.tf_to_motif[x].values()))
        # Keep the highest IC amongst the duplicates
        to_drop = duplicated['Original Name'][duplicated.sort_values(
            'IC').duplicated('Human Name', keep='last')]

        # Exlude the IDs which are already present among the human ones
        human_motif_names = [i.name for i in self.motifs if 
            '9606' in i.species]
        corr_df['Duplicate'] = corr_df['Human Name'].isin(human_motif_names)

        # Perform the final filtering - discard all duplicates and TFs without 
        # homologs
        corr_df_final = corr_df[(corr_df['Duplicate'] == False) & 
            (corr_df['Human Name'] != '') & 
            (corr_df['Original Name'].isin(to_drop) == False)]

        # The mapping of original to human names and the matrix IDs to be kept
        animal_to_human = {i: j for i,j in zip(corr_df_final['Original Name'], 
            corr_df_final['Human Name'])}
        print ()
        print ('Final number of IDs which will be replaced by human homologs:', 
               len(animal_to_human))
        vert_list_nohd = [i.matrix_id for i in self.motifs if 
            (i.name in human_motif_names or i.name in animal_to_human)]
        print ('Final number of total matrix IDs:', len(vert_list_nohd))

        self.matrix_ids = vert_list_nohd
        self.tf_names = human_motif_names + list(animal_to_human.keys())
        self.animal_to_human = animal_to_human


    def filter_matches(
        self, 
        promoter_file=None, 
        bigbed_file=None,
        n_processes=1
    ):
        
        # TODO: Implement defaults and downloads
        if promoter_file is None:
            promoter_file = self.retrieve_file('promoter', self.temp_folder)
            if promoter_file is None:
                print ('Unable to find or retrieve the promoter file, exiting')
                return
        if bigbed_file is None:
            bigbed_file = self.retrieve_file('jaspar_bigbed', self.temp_folder)
            if bigbed_file is None:
                print ('Unable to find or retrieve the JASPAR bigbed file, exiting')
                return

        print ('Loading the promoter bed file...')
        print ()
        df_full = bioframe.read_table(promoter_file, schema='bed')
        df_full['name'] = df_full['name'].apply(lambda x: x.split('.')[0])
        df_full.drop(columns=['score', 'strand'], inplace=True)
        df_full.set_index('name', inplace=True)

        chr_list = ['chr{}'.format(i) for i in [j for j in range(1,23)] + 
            ['M', 'X', 'Y']]
        results_list = []
        p = Pool(n_processes)

        print ('Iterating over the chromosomes...')
        start_time = time.time()
        for chrom in chr_list:
            st_chr = time.time()
            df_chrom = df_full[df_full['chrom'] == chrom]
            if len(df_chrom) == 0:
                suffix = 'no transcripts'
            elif len(df_chrom) == 1:
                suffix = '1 transcript'
            else:
                suffix = f'{len(df_chrom)} transcripts'
            print (f'Chromosome {chrom[3:]} with ' + suffix)
            chunk_size = ceil(sqrt(len(df_chrom) / n_processes))
            chunk_divisions = [i for i in range(0, len(df_chrom), chunk_size)]
            input_tuples = [(bigbed_file, df_chrom, self.tf_names, chrom, i, i+chunk_size)
                for i in chunk_divisions]
            result = p.map_async(filter_edges_helper, input_tuples, chunksize=n_processes)
            edges_chrom_list = result.get()
            results_list += edges_chrom_list
            elapsed_chr = time.time() - st_chr
            print (f'Done in: {elapsed_chr // 60:n} m {elapsed_chr % 60:.2f} s')

        elapsed = time.time() - start_time
        print ()
        print (f'Total time: {elapsed // 60:n} m {elapsed % 60:.2f} s')

        self.all_edges = pd.concat(results_list, ignore_index=True)