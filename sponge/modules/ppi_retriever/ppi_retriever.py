### Imports ###
import requests

import numpy as np
import pandas as pd

from io import BytesIO
from typing import Iterable

from sponge.config_manager import ConfigManager
from sponge.modules.motif_selector import ProteinIDMapper
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class PPIRetriever:
    # Functions
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):
        
        self.core_config = core_config
        self.ppi_url = core_config['url']['ppi']
        self.protein_url = core_config['url']['protein']
        self.version_logger = version_logger
        self.physical_only = not user_config.is_false(['ppi', 'physical_only'])


    def retrieve_ppi(
        self,
        tf_names: Iterable[str],
    ):
        """
        Retrieves the protein-protein interaction data from the STRING
        database for the previously identified transcription factors.
        Stores the resulting network internally.
        """

        print ()
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
        self.version_logger.write_retrieved('string_ppi',
            version_df['string_version'].loc[0])

        if len(ids_to_check) > 0:
            # Retrieve UniProt identifiers for the genes with differing names
            print ('Checking the conflicts in the UniProt database...')
            mapper = ProteinIDMapper(self.core_config)
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
        query_string_filt = '%0d'.join(matching_ids)

        print ('Retrieving the network from STRING...')
        network_str = (f'{self.ppi_url}network?'
            f'identifiers={query_string_filt}&species=9606')
        if self.physical_only:
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

        print ()
        print ('Final number of TFs in the PPI network: '
            f'{len(set(ppi_df["tf1"]).union(set(ppi_df["tf2"])))}')
        print (f'Final number of edges: {len(ppi_df)}')

        self.ppi_frame = ppi_df