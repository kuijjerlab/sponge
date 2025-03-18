### Imports ###
import requests

import pandas as pd

from sponge.config_manager import ConfigManager
from sponge.modules.utils import adjust_gene_name
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class HomologyRetriever:
    # Functions
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):

        self.unique_motifs = user_config['motif']['unique_motifs']
        self.ncbi_url = core_config['url']['homology']
        self.version_logger = version_logger


    def find_homologs(
        self,
        motifs,
        tf_to_motif,
        mapping,
    ) -> None:

        # TODO: Replace mentions of human to make later expansion easier
        # Get the non-human motifs
        non_human_motifs = [i for i in motifs if '9606' not in i.species]
        print ()
        print ('Non-human motifs:', len(non_human_motifs))

        # Get the human homologs from NCBI
        print ('Retrieving homologs from NCBI...')
        homologs = {}
        suffix = '/gene/id/{gene_id}/orthologs'
        for acc,gene_id in mapping[['from', 'to']].values:
            r = requests.get(self.ncbi_url + suffix.format(gene_id=gene_id),
                params=dict(taxon_filter=9606))
            r.raise_for_status()
            table = r.json()
            if 'reports' in table:
                gene = table['reports'][0]['gene']
                homologs[acc] = [gene['symbol'], gene['gene_id']]
        # Record the version of NCBI services
        version_r = requests.get(self.ncbi_url + '/version')
        version = version_r.json()['version']
        self.version_logger.write_retrieved('NCBI', version)

        # Get the non-human motif names
        non_human_motif_names = [i.name for i in non_human_motifs]
        # Compare against NCBI homologs
        found_names = [adjust_gene_name(i.name) for i in non_human_motifs
            if not False in [acc in homologs for acc in i.acc]]
        # Find the missing ones
        missing = (set([adjust_gene_name(i) for i in non_human_motif_names]) -
            set(found_names))
        print ()
        print ('TFs for which no homolog was found:')
        for motif in non_human_motifs:
            if motif.name in missing:
                print (motif.name, *motif.acc)

        # Create a DataFrame of corresponding names
        corr_names = {
            motif.name: '::'.join([homologs[acc][0] for acc in motif.acc])
            for motif in non_human_motifs
            if motif.name in found_names
        }
        corr_df = pd.DataFrame(non_human_motif_names,
            columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(
            adjust_gene_name)
        corr_df['Human Name'] = corr_df['Original Name'].apply(corr_names.get)

        if self.unique_motifs:
            # Find duplicates
            duplicated = corr_df[corr_df['Human Name'].duplicated(keep=False) &
                ~corr_df['Human Name'].isna()].copy()
            to_print = duplicated.groupby('Human Name'
                )['Original Name'].unique().apply(lambda x: ' '.join(x))
            print ()
            print ('Duplicate names:')
            for i in to_print.index:
                print (f'{i}:', to_print.loc[i])

            # Calculate the information content for duplicates
            duplicated['IC'] = duplicated['Original Name'].apply(lambda x:
                max(tf_to_motif[x].values()))
            # Keep the highest IC amongst the duplicates
            to_drop = duplicated['Original Name'][duplicated.sort_values(
                'IC').duplicated('Human Name', keep='last')]
        else:
            to_drop = []

        # Exclude the IDs which are already present among the human ones
        human_motif_names = [i.name for i in motifs if '9606' in i.species]
        corr_df['Duplicate'] = corr_df['Human Name'].isin(human_motif_names)

        # Perform the final filtering - discard all duplicates and TFs without
        # homologs
        corr_df_final = corr_df[(corr_df['Duplicate'] == False) &
            (~corr_df['Human Name'].isna()) &
            (corr_df['Original Name'].isin(to_drop) == False)]

        # The mapping of original to human names and the matrix IDs to be kept
        animal_to_human = {animal_name: human_name for animal_name, human_name
            in zip(corr_df_final['Original Name'],
                corr_df_final['Human Name'])}
        print ()
        print ('Final number of IDs which will be replaced by homologs:',
               len(animal_to_human))
        # Doing it this way ensures the ordering matches
        matrix_ids = [motif.matrix_id for motif in motifs if
            (motif.name in human_motif_names or motif.name in animal_to_human)]
        tf_names = [motif.name for motif in motifs if
            (motif.name in human_motif_names or motif.name in animal_to_human)]
        print ('Final number of all matrix IDs:', len(matrix_ids))

        self.animal_to_human = animal_to_human
        self.matrix_ids = matrix_ids
        self.tf_names = tf_names