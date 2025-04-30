### Imports ###
import requests
import time

import pandas as pd

from typing import Iterable, Union

### Class definition ###
class ProteinIDMapper:
    # Methods
    def __init__(
        self,
        mapping_url: str,
    ):
        """
        Class that maps protein IDs using the UniProt services.

        Parameters
        ----------
        mapping_url : str
            URL of the UniProt ID mapping service
        """

        self.mapping_url = mapping_url


    def get_uniprot_mapping(
        self,
        from_db: str,
        to_db: str,
        ids: Union[str, Iterable[str]],
        **kwargs,
    ) -> pd.DataFrame:
        """
        Attempts to get a mapping for the given IDs from UniProt. Can be
        provided with extra keyword arguments which are then added to
        the request.

        Parameters
        ----------
        from_db : str
            Name of the database to match from
        to_db : str
            Name of the database to match to
        ids : Union[str, Iterable[str]]
            Single ID or Iterable of IDs to match

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame containing the mapping

        Raises
        ------
        requests.exceptions.HTTPError
            Reproduction of an error message from UniProt if no job ID
            was retrieved, typically pointing to an issue with the query
        """

        # Guard against empty request
        if len(ids) == 0:
            return pd.DataFrame(columns=['from', 'to'])
        # The basic form of the request
        data = {'ids': ids, 'from': from_db, 'to': to_db}
        # Potential additional arguments
        data.update(kwargs)
        # Post the request and register the reply
        uniprot_request = requests.post(self.mapping_url + 'run', data)
        uniprot_reply = uniprot_request.json()
        if 'jobId' in uniprot_reply:
            job_id = uniprot_reply['jobId']
        else:
            # No job ID was assigned - probably an issue with the query
            raise requests.exceptions.HTTPError(uniprot_reply['messages'][0])

        MAX_ITERATIONS = 40
        for _ in range(MAX_ITERATIONS):
            # Loop until the results are available or max iterations exceeded
            uniprot_status = requests.get(self.mapping_url +
                f'status/{job_id}')
            if 'results' in uniprot_status.json():
                break
            # Try every half a second
            time.sleep(0.5)
        if 'results' not in uniprot_status.json():
            # Unable to retrieve the results within the given time
            print ('No results have been retrieved in the given time')
            return pd.DataFrame(columns=['from', 'to'])

        # Retrieve the results
        uniprot_results = requests.get(self.mapping_url + f'stream/{job_id}')

        # Convert the results to a pandas DataFrame
        results_df = pd.DataFrame(uniprot_results.json()['results'])
        results_df.drop_duplicates(subset='from', inplace=True)

        return results_df