### Imports ###
import requests
import time
import os

from biomart import BiomartServer

from typing import Optional, Union, Iterable, Tuple

from io import BytesIO

from collections import defaultdict

from tqdm import tqdm

from sponge.config import *

### Functions ###
def prompt_to_confirm(
    question: str
) -> bool:
    
    key = None
    positive = ['y', 'yes', 'hell yeah']
    negative = ['n', 'no', 'nope']
    while key is None or key.lower() not in positive + negative:
        if key is not None:
            print (f'Input not recognised: {key}')
        print (f'{question} Y/N', flush=True)
        key = input()
        print (key)
    print ()

    return key.lower() in positive


def description_to_path(
    description: str,
    temp_folder: Union[str, bytes, os.PathLike]
) -> Optional[str]:

    if description not in FILE_DF.index:
        print (f'File description not recognised: {description}')
        return None
    file_name = FILE_DF.loc[description, 'name']
    file_path = os.path.join(temp_folder, file_name)

    return file_path


def check_file_exists(
    description: str
) -> bool:
    
    return os.path.exists(description_to_path(description))


def load_promoters_from_biomart(
    file_path: Union[str, bytes, os.PathLike],
    filter_basic: bool = True,
    chromosomes: Iterable[str] = [str(i) for i in range(1,23)] + 
        ['MT', 'X', 'Y'],
    tss_offset: Tuple[int, int] = (-750, 250),
    keep_ensembl: bool = True
) -> dict:

    answer = {}
    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_transcript_id', 'transcript_gencode_basic', 
        'chromosome_name', 'transcription_start_site', 'strand']
    if keep_ensembl:
        attributes += ['ensembl_gene_id', 'external_gene_name', 
            'gene_biotype']
    print ('Retrieving response to query...')
    response = ensembl.search({'attributes': attributes}, header=1)
    buffer = download_with_progress(response)
    answer['version'] = ensembl.display_name
    dtype_dict = defaultdict(lambda: str)
    dtype_dict['Transcription start site (TSS)'] = int
    dtype_dict['Strand'] = int
    df = pd.read_csv(buffer, sep='\t', dtype=dtype_dict)
    print ('Filtering and modifying dataframe...')
    if filter_basic:
        df = df[df['GENCODE basic annotation'] == 'GENCODE basic'].copy()
        df.drop(columns='GENCODE basic annotation', inplace=True)
    if chromosomes is not None:
        df = df[df['Chromosome/scaffold name'].isin(chromosomes)]
    df['Chromosome'] = df['Chromosome/scaffold name'].apply(lambda x: 
        'chrM' if x == 'MT' else f'chr{x}')
    df['Strand'] = df['Strand'].apply(lambda x: '+' if x > 0 else '-')
    df['Start'] = df.apply(lambda row: 
        row['Transcription start site (TSS)'] + tss_offset[0] 
        if row['Strand'] == '+' 
        else row['Transcription start site (TSS)'] - tss_offset[1], 
        axis=1)
    df['End'] = df['Start'] + (tss_offset[1] - tss_offset[0])
    df['Score'] = 0
    df.sort_values(['Chromosome', 'Start'], inplace=True)
    columns = ['Chromosome', 'Start', 'End', 'Transcript stable ID', 
        'Score', 'Strand']
    print (f'Saving data to {file_path}...')
    df[columns].to_csv(file_path, sep='\t', header=False, index=False)
    print ()
    if keep_ensembl:
        answer['ensembl'] = df[['Gene stable ID', 'Transcript stable ID', 
            'Gene name', 'Gene type']]

    return answer


def load_ensembl_from_biomart(
    file_path: Union[str, bytes, os.PathLike]
) -> None:
    
    answer = {}
    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 
        'external_gene_name', 'gene_biotype']
    print ('Retrieving response to query...')
    response = ensembl.search({'attributes': attributes}, header=1)
    buffer = download_with_progress(response)
    answer['version'] = ensembl.display_name
    df = pd.read_csv(buffer, sep='\t')     
    df.to_csv(file_path, sep='\t', index=False)
    answer['ensembl'] = df

    return answer


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
    MAX_ITERATIONS = 40
    for _ in range(MAX_ITERATIONS):
        # Loop until the results are available
        uniprot_status = requests.get(MAPPING_URL + f'status/{job_id}')
        if 'results' in uniprot_status.json():
            break
        # Try every half a second
        time.sleep(0.5)
    if 'results' not in uniprot_status.json():
        # Unable to retrieve the results within the given time
        print ('No results have been retrieved in the given time')
        return pd.DataFrame()
    # Retrieve the results
    uniprot_results = requests.get(MAPPING_URL + f'stream/{job_id}')
    # Convert the results to a pandas DataFrame
    results_df = pd.DataFrame(uniprot_results.json()['results'])
    results_df.drop_duplicates(subset='from', inplace=True)

    return results_df