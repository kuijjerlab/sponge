### Imports ###
import bioframe
import requests
import time
import os

from typing import Optional, Union, Iterable, Tuple

from io import BytesIO

from tqdm import tqdm

from sponge.config import *

# TODO: Move file related stuff into this file from Sponge if possible

### Functions ###
def download_with_progress(
    url: Union[str, requests.models.Response],
    file_path: Optional[Union[str, bytes, os.PathLike]] = None,
    desc: str = 'response'
) -> Optional[BytesIO]:
    """
    Downloads from a given URL or retrieves a response to a given 
    request while providing a progress bar.

    Parameters
    ----------
    url : Union[str, requests.models.Response]
        The URL or response to be processed
    file_path : Optional[Union[str, bytes, os.PathLike]], optional
        The file path for saving or None to save into a BytesIO object,
        by default None
    desc : str, optional
        The description to show, by default 'response'

    Returns
    -------
    Optional[BytesIO]
        The BytesIO object containing the data or None if file_path was 
        not set to None
    """
    
    # Determine the type of request
    if type(url) == str:
        request = requests.get(url, stream=True)
    else:
        request = url
    total = int(request.headers.get('content-length', 0))   
    # Determine whether to save data to a file or object       
    if file_path is None:
        stream = BytesIO()
    else:
        stream = open(file_path, 'wb')
        desc = file_path
    # Download with a progress bar using tqdm
    with tqdm(desc=desc, total=total, unit='iB', unit_scale=True,
        unit_divisor=1024) as bar:
        for data in request.iter_content(chunk_size=1024):
            size = stream.write(data)
            bar.update(size)
    
    if file_path is None:
        return BytesIO(stream.getvalue())


def get_uniprot_mapping(
    from_db: str,
    to_db: str,
    ids: Union[str, Iterable[str]],
    **kwargs
) -> pd.DataFrame:
    """
    Attempts to get a mapping for the given IDs from Uniprot. Can be 
    provided with extra keyword arguments which are then added to the 
    request.

    Parameters
    ----------
    from_db : str
        The name of the database to match from
    to_db : str
        The name of the database to match to
    ids : Union[str, Iterable[str]]
        A single ID or an Iterable of IDs to match

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the mapping

    Raises
    ------
    RuntimeError
        Reproduction of an error message from UniProt if no job ID
        was retrieved, typically pointing to an issue with the query
    """

    # TODO: Uniprot has a limit on the number of ids that can be checked at one
    # time, test if this is affected  
    # The basic form of the request
    data = {'ids': ids, 'from': from_db, 'to': to_db}
    # Potential additional arguments
    data.update(kwargs)
    # Post the request and register the reply
    uniprot_request = requests.post(MAPPING_URL + 'run', data)
    uniprot_reply = uniprot_request.json()
    if 'jobId' in uniprot_reply:
        job_id = uniprot_reply['jobId']
    else:
        # No job ID was assigned - probably an issue with the query
        print ('Unable to retrieve a job ID from UniProt')
        print ('The following reply was received instead:')
        for message in uniprot_reply['messages']:
            print (message)
        raise RuntimeError()
    # TODO: Potentially implement max_iterations
    while True:
        # Loop until the results are available
        uniprot_status = requests.get(MAPPING_URL + f'status/{job_id}')
        if 'results' in uniprot_status.json():
            break
        # Try every half a second
        time.sleep(0.5)
    # Retrieve the results
    uniprot_results = requests.get(MAPPING_URL + f'stream/{job_id}')
    # Convert the results to a pandas DataFrame
    results_df = pd.DataFrame(uniprot_results.json()['results'])
    results_df.drop_duplicates(subset='from', inplace=True)

    return results_df


def filter_edges(
    bb_ref: Union[str, bytes, os.PathLike], 
    bed_df: pd.DataFrame, 
    motif_list: Iterable[str], 
    chrom: str, 
    start_ind: int, 
    final_ind: int, 
    score_threshold: float = 400
) -> pd.DataFrame:
    """
    Filters possible binding site matches for the provided transcription
    factors from a bigbed file. This is done for a number of continuous
    regions of a single chromosome subject to a score threshold.

    Parameters
    ----------
    bb_ref : Union[str, bytes, os.PathLike]
        The path to a bigbed file that stores all possible matches
    bed_df : pd.DataFrame
        A pandas DataFrame containing the regions of interest in the
        chromosome, typically promoters
    motif_list : Iterable[str]
        An Iterable containing the names of transcription factors of
        interest
    chrom : str
        The name of the chromosome of interest
    start_ind : int
        The starting index of the region DataFrame (bed_df)
    final_ind : int
        The final index of the region DataFrame (bed_df)
    score_threshold : float, optional
        The score required for selection, by default 400

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered edges from the 
        regions of interest
    """
    
    df = pd.DataFrame()
    for transcript in bed_df.index[start_ind:final_ind]:
        # Retrieve the start and end points
        start,end = bed_df.loc[transcript][['start', 'end']]
        # Load all matches in that region from the bigbed file
        motifs = bioframe.read_bigbed(bb_ref, chrom, start=start, end=end)
        # Filter only the transcription factors in the list
        max_scores = motifs[motifs['TFName'].isin(motif_list)].groupby(
            'TFName')[['score', 'name']].max()
        # Filter only high enough scores
        max_scores = max_scores[max_scores['score'] >= score_threshold]
        max_scores.reset_index(inplace=True)
        # Add the transcript (region) name for identification
        max_scores['transcript'] = transcript
        # Append the results to the dataframe
        df = pd.concat((df, max_scores), ignore_index=True, copy=False)

    return df


def filter_edges_helper(
    input_tuple: Tuple[str, pd.DataFrame, Iterable[str], str, int, int, float]
) -> pd.DataFrame:
    """
    Serves as a wrapper around the filter_edges function that only
    takes a single input, allowing it to be processed by map_async from
    the multiprocessing module.

    Parameters
    ----------
    input_tuple : Tuple[str, pd.DataFrame, Iterable[str], str, int, int, 
        float]
        A tuple of all inputs to the filter_edges function, for more
        details refer to its docstring

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered edges corresponding
        to the given input
    """

    return filter_edges(*input_tuple)