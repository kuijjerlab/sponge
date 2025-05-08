### Imports ###
import requests

import numpy as np
import pandas as pd

from io import BytesIO
from typing import Iterable

from sponge.modules.protein_id_mapper import ProteinIDMapper

### Class definition ###
class PPIRetriever:
    # Methods
    def __init__(
        self,
        ppi_url: str,
        protein_url: str,
    ):
        """
        Class that retrieves PPI data from the STRING database.

        Parameters
        ----------
        ppi_url : str
            URL of the STRING database API
        protein_url : str
            URL of the UniProt ID mapping service
        """

        self.ppi_url = ppi_url
        self.protein_url = protein_url
        # To be overwritten once retrieved
        self.ppi_frame = None
        # Overwritten by registering with a VersionLogger instance
        self.version_logger = None


    def retrieve_ppi(
        self,
        tf_names: Iterable[str],
        score_threshold: int = 400,
        physical_only: bool = True,
    ):
        """
        Retrieves the protein-protein interaction data from the STRING
        database for the provided proteins. Stores the resulting network
        internally.

        Parameters
        ----------
        tf_names : Iterable[str]
            Names of the proteins
        score_threshold : float, optional
            Minimal interaction score for it to be included in the
            network, by default 400
        physical_only : bool, optional
            Whether to consider only physical interactions,
            by default True
        """

        print ('\n--- Retrieving protein-protein interaction data ---')

        if len(tf_names) <= 1:
            print ('Not enough TFs were provided to build a PPI network.')
            self.ppi_frame = pd.DataFrame(columns=['tf1', 'tf2', 'score'])
            return

        print ('Retrieving mapping from STRING...')
        query_string = '%0d'.join(tf_names)
        mapping_request = requests.get(f'{self.ppi_url}get_string_ids?'
            f'identifiers={query_string}&species=9606')
        mapping_df = pd.read_csv(BytesIO(mapping_request.content), sep='\t')
        mapping_df['queryName'] = mapping_df['queryIndex'].apply(
            lambda i: tf_names[i])
        # Check where the preferred name doesn't match the query
        diff_df = mapping_df[mapping_df['queryName'] !=
            mapping_df['preferredName']]
        ids_to_check = np.concatenate((diff_df['queryName'],
            diff_df['preferredName']))
        matching_ids = list(mapping_df[mapping_df['queryName'] ==
            mapping_df['preferredName']]['preferredName'])
        # Log the STRING version in the fingerprint
        version_request = requests.get(f'{self.ppi_url}version')
        version_df = pd.read_csv(BytesIO(version_request.content), sep='\t',
            dtype=str)
        self.write_retrieved('string_ppi', version_df['string_version'].loc[0])

        if len(ids_to_check) > 0:
            # Retrieve UniProt identifiers for the genes with differing names
            print ('Checking the conflicts in the UniProt database...')
            mapper = ProteinIDMapper(self.protein_url)
            uniprot_df = mapper.get_uniprot_mapping('Gene_Name', 'UniProtKB',
                ids_to_check).set_index('from')
            p_to_q = {p: q for q,p in zip(diff_df['queryName'],
                diff_df['preferredName'])}
            # Keep the conflicts where there is a match or where one or both
            # of the names doesn't find an identifier
            for p,q in p_to_q.items():
                if (p not in uniprot_df.index or q not in uniprot_df.index
                    or uniprot_df.loc[p, 'to'] == uniprot_df.loc[q, 'to']):
                    matching_ids.append(p)

        if len(matching_ids) <= 1:
            print ('Not enough IDs were matched to build a PPI network.')
            self.ppi_frame = pd.DataFrame(columns=['tf1', 'tf2', 'score'])
            return

        print ('Retrieving the network from STRING...')
        query_string_filt = '%0d'.join(matching_ids)
        network_str = (f'{self.ppi_url}network?'
            f'identifiers={query_string_filt}&species=9606&'
            f'required_score={score_threshold}')
        if physical_only:
            network_str += '&network_type=physical'
        request = requests.get(network_str)
        ppi_df = pd.read_csv(BytesIO(request.content), sep='\t')

        print ('Processing the results...')
        ppi_df.drop(['stringId_A', 'stringId_B', 'ncbiTaxonId', 'nscore',
            'fscore', 'pscore', 'ascore', 'escore', 'dscore', 'tscore'],
            axis=1, inplace=True)
        ppi_df.rename(columns={'preferredName_A': 'tf1',
            'preferredName_B': 'tf2'}, inplace=True)
        if len(ids_to_check) > 0:
            # Replace with names that have been queried (as used by JASPAR)
            ppi_df['tf1'].replace(p_to_q, inplace=True)
            ppi_df['tf2'].replace(p_to_q, inplace=True)
        ppi_df.sort_values(by=['tf1', 'tf2'], inplace=True)

        print ('\nFinal number of TFs in the PPI network: '
            f'{len(set(ppi_df["tf1"]).union(set(ppi_df["tf2"])))}')
        print (f'Final number of edges: {len(ppi_df)}')

        self.ppi_frame = ppi_df


    def get_ppi_frame(
        self,
    ) -> pd.DataFrame:
        """
        Get the pandas DataFrame of the PPI network, can be None
        if it wasn't retrieved yet.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the PPI network
        """

        return self.ppi_frame

    # Placeholder functions to be replaced with VersionLogger if any
    def write_retrieved(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass