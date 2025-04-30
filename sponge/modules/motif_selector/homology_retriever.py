### Imports ###
import requests

import pandas as pd

from Bio.motifs.jaspar import Motif
from typing import List, Mapping

from sponge.modules.protein_id_mapper import ProteinIDMapper
from sponge.modules.utils import adjust_gene_name

### Class definition ###
class HomologyRetriever:
    # Methods
    def __init__(
        self,
        unique_motifs: bool,
        mapping_url: str,
        ncbi_url: str,
    ):
        """
        Class that retrieves TF homologs from NCBI.

        Parameters
        ----------
        unique_motifs : bool
            Whether a single motif is enforced per TF
        mapping_url : str
            URL to the ID mapping service of UniProt
        ncbi_url : str
            URL to the NCBI datasets API
        """

        self.unique_motifs = unique_motifs
        self.mapping_url = mapping_url
        self.ncbi_url = ncbi_url
        # To be overwritten once retrieved
        self.homolog_mapping = None
        self.matrix_ids = None
        self.tf_names = None
        # Overwritten by registering with a VersionLogger instance
        self.version_logger = None


    def find_homologs(
        self,
        motifs: List[Motif],
        tf_to_motif: Mapping[str, dict],
    ) -> None:
        """
        Attempts to find the homologs of the provided motifs that do not
        match the species of interest.

        Parameters
        ----------
        motifs : List[Motif]
            List of selected TF motifs
        tf_to_motif : Mapping[str, dict]
            Mapping of TF names to their motif IDs and information
            contents
        """

        print ('\n--- Searching for matching homologs ---')

        # Get the non-species motifs
        xeno_motifs = [motif for motif in motifs
            if '9606' not in motif.species]
        print ('Motifs from other species:', len(xeno_motifs))

        # Retrieve mapping of UniProt to GeneID
        all_ids = set()
        for motif in xeno_motifs:
            for id in motif.acc:
                all_ids.add(id)
        mapper = ProteinIDMapper(self.mapping_url)
        mapping = mapper.get_uniprot_mapping('UniProtKB_AC-ID', 'GeneID',
            list(all_ids))

        # Get the homologs from NCBI
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
        self.write_retrieved('ncbi_homologs', version)

        # Get the non-species motif names
        xeno_motif_names = [motif.name for motif in xeno_motifs]
        # Compare against NCBI homologs
        found_names = [adjust_gene_name(motif.name) for motif in xeno_motifs
            if not False in [acc in homologs for acc in motif.acc]]
        # Find the missing ones
        missing = (set([adjust_gene_name(name) for name in
            xeno_motif_names]) - set(found_names))
        print ('\nTFs for which no homolog was found:')
        for motif in xeno_motifs:
            if motif.name in missing:
                print (motif.name, *motif.acc)

        # Create a DataFrame of corresponding names
        corr_names = {
            motif.name: '::'.join([homologs[acc][0] for acc in motif.acc])
            for motif in xeno_motifs
            if not False in [acc in homologs for acc in motif.acc]
        }
        corr_df = pd.DataFrame(xeno_motif_names, columns=['Original Name'])
        corr_df['Adjusted Name'] = corr_df['Original Name'].apply(
            adjust_gene_name)
        corr_df['Species Name'] = corr_df['Original Name'].apply(
            corr_names.get)

        if self.unique_motifs:
            # Find duplicates
            duplicated = corr_df[corr_df['Species Name'].duplicated(keep=False
                ) & ~corr_df['Species Name'].isna()].copy()
            to_print = duplicated.groupby('Species Name'
                )['Original Name'].unique().apply(lambda x: ' '.join(x))
            print ('\nDuplicate names:')
            for i in to_print.index:
                print (f'{i}:', to_print.loc[i])

            # Calculate the information content for duplicates
            duplicated['IC'] = duplicated['Original Name'].apply(lambda x:
                max(tf_to_motif[x].values()))
            # Keep the highest IC amongst the duplicates
            to_drop = duplicated['Original Name'][duplicated.sort_values(
                'IC').duplicated('Species Name', keep='last')]
        else:
            to_drop = []

        # Exclude the IDs which are already present among the in-species ones
        species_motif_names = [motif.name for motif in motifs
            if '9606' in motif.species]
        corr_df['Duplicate'] = corr_df['Species Name'].isin(
            species_motif_names)

        # Perform the final filtering - discard all duplicates and TFs without
        # homologs
        corr_df_final = corr_df[(corr_df['Duplicate'] == False) &
            (~corr_df['Species Name'].isna()) &
            (corr_df['Original Name'].isin(to_drop) == False)]

        # The mapping of out-species to in-species names
        # and the matrix IDs to be kept
        homolog_mapping = {animal_name: species_name
            for animal_name,species_name
            in zip(corr_df_final['Original Name'],
                corr_df_final['Species Name'])}
        print ('\nFinal number of IDs which will be replaced by homologs:',
               len(homolog_mapping))
        # Doing it this way ensures the ordering matches
        matrix_ids = [motif.matrix_id for motif in motifs if
            (motif.name in species_motif_names or
            motif.name in homolog_mapping)]
        tf_names = [motif.name for motif in motifs if
            (motif.name in species_motif_names or
            motif.name in homolog_mapping)]
        print ('Final number of all matrix IDs:', len(matrix_ids))

        self.homolog_mapping = homolog_mapping
        self.matrix_ids = matrix_ids
        self.tf_names = tf_names


    def get_homolog_mapping(
        self,
    ) -> Mapping[str, str]:
        """
        Get the mapping of TF homolog names to the names in the species
        of interest, can be None if it wasn't retrieved.

        Returns
        -------
        Mapping[str, str]
            Mapping of TF homolog names
        """

        return self.homolog_mapping


    def get_matrix_ids(
        self,
    ) -> List[str]:
        """
        Get the list of TF matrix IDs, can be None if it wasn't
        retrieved.

        Returns
        -------
        List[str]
            List of TF matrix IDs
        """

        return self.matrix_ids


    def get_tf_names(
        self,
    ) -> List[str]:
        """
        Get the list of TF names, can be None if it wasn't retrieved.

        Returns
        -------
        List[str]
            List of TF names
        """

        return self.tf_names

    # Placeholder functions to be replaced with VersionLogger if any
    def write_retrieved(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass