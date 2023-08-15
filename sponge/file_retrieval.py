### Imports ###
import requests
import time
import os
import gzip

from biomart import BiomartServer

from typing import Optional, Union, Iterable, Tuple, Dict

from io import BytesIO

from collections import defaultdict

from tqdm import tqdm

from sponge.config import *

PATH = Union[str, os.PathLike]
FILE_LIKE = Union[str, bytes, os.PathLike]

### Functions ###
def prompt_to_confirm(
    question: str
) -> bool:
    """
    Asks the user to confirm the choice interactively, using the
    provided question.

    Parameters
    ----------
    question : str
        The yes or no question to answer

    Returns
    -------
    bool
        The response provided by user
    """
    
    key = None
    # Accepted replies
    positive = ['y', 'yes', 'hell yeah']
    negative = ['n', 'no', 'nope']
    while key is None or key.lower() not in positive + negative:
        if key is not None:
            # This means there was a reply, but not an accepted one
            print (f'Input not recognised: {key}')
        print (f'{question} Y/N', flush=True)
        key = input()
        print (key)
    print ()

    return key.lower() in positive


def description_to_path(
    description: str,
    temp_folder: PATH
) -> Optional[str]:
    """
    Converts the description of a file to its expected path.

    Parameters
    ----------
    description : str
        The description of a file
    temp_folder : PATH
        The path to the temp folder where files are located

    Returns
    -------
    Optional[str]
        The expected path to the file, or None if description is not
        recognised
    """

    # All the valid descriptions are indices of the file DataFrame
    if description not in FILE_DF.index:
        print (f'File description not recognised: {description}')
        return None
    # Join the folder and the expected name
    file_name = FILE_DF.loc[description, 'name']
    file_path = os.path.join(temp_folder, file_name)

    return file_path


def check_file_exists(
    description: str,
    temp_folder: PATH
) -> bool:
    """
    Checks if the file corresponding to the description exists.

    Parameters
    ----------
    description : str
        The description of a file
    temp_folder : PATH
        The path to the temp folder where files are located

    Returns
    -------
    bool
        Whether the file corresponding to the description exists
    """
    
    return os.path.exists(description_to_path(description, temp_folder))


def load_promoters_from_biomart(
    file_path: FILE_LIKE,
    filter_basic: bool = True,
    chromosomes: Optional[Iterable[str]] = 
        [str(i) for i in range(1,23)] + ['MT', 'X', 'Y'],
    chromosome_mapping: pd.Series = DEFAULT_MAPPING,
    tss_offset: Tuple[int, int] = (-750, 250),
    keep_ensembl: bool = True
) -> Dict[str, Union[str, pd.DataFrame]]:
    """
    Generates the promoter file from the data retrieved from the Ensembl
    BioMart server. Optionally also keeps a subset of the data as a 
    DataFrame for downstream use. 

    Parameters
    ----------
    file_path : FILE_LIKE
        The path to where the resulting file should be saved
    filter_basic : bool, optional
        Whether to filter for only the GENCODE basic transcripts, 
        by default True
    chromosomes : Optional[Iterable[str]], optional
        An Iterable of chromosomes to be considered or None to consider 
        all, by default [str(i) for i in range(1,23)] + ['MT', 'X', 'Y']
    chromosome_mapping : pd.Series, optional
        The mapping of Ensembl chromosome names to the UCSC ones, by 
        default a simple mapping of only the main chromosomes and MT
    tss_offset : Tuple[int, int], optional
        The offset from the transcription start site to define the 
        promoter region, by default (-750, 250)
    keep_ensembl : bool, optional
        Whether to return the Ensembl DataFrame with a subset of the 
        data (gene and transcript IDs, gene name, gene type), 
        by default True

    Returns
    -------
    Dict[str, Union[str, pd.DataFrame]]
        A dictionary containing the version of the database used and
        optionally the Ensembl DataFrame
    """

    answer = {}

    # Select the right dataset from BioMart
    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    # Attributes to retrieve
    attributes = ['ensembl_transcript_id', 'transcript_gencode_basic', 
        'chromosome_name', 'transcription_start_site', 'strand']
    # Extra attributes that matter only for the Ensembl DataFrame
    if keep_ensembl:
        attributes += ['ensembl_gene_id', 'external_gene_name', 
            'gene_biotype']
    print ('Retrieving response to query...')
    # Submit and retrieve the response
    response = ensembl.search({'attributes': attributes}, header=1)
    buffer = download_with_progress(response)

    # Save the database version into the dictionary
    answer['version'] = ensembl.display_name
    # Dictionary of types for conversion from the response, default strings
    dtype_dict = defaultdict(lambda: str)
    # Change the types that are not strings but integers
    dtype_dict['Transcription start site (TSS)'] = int
    dtype_dict['Strand'] = int
    # Convert the response into a DataFrame
    df = pd.read_csv(buffer, sep='\t', dtype=dtype_dict)
    
    print ('Filtering and modifying dataframe...')
    if filter_basic:
        # Filter only for GENCODE basic
        df = df[df['GENCODE basic annotation'] == 'GENCODE basic'].copy()
        df.drop(columns='GENCODE basic annotation', inplace=True)
    if chromosomes is not None:
        # Filter only for selected chromosomes
        df = df[df['Chromosome/scaffold name'].isin(chromosomes)]
    # Convert chromosome names to match with other inputs
    df['Chromosome'] = df['Chromosome/scaffold name'].apply(lambda x: 
        chromosome_mapping[x])
    # Convert strand to +/-
    df['Strand'] = df['Strand'].apply(lambda x: '+' if x > 0 else '-')
    # Calculate the start based on the given offset from TSS
    # The calculation is strand dependent
    df['Start'] = df.apply(lambda row:
        row['Transcription start site (TSS)'] + tss_offset[0]
        if row['Strand'] == '+'
        else row['Transcription start site (TSS)'] - tss_offset[1], 
        axis=1)
    # End is always greater than start, this way it is strand independent
    df['End'] = df['Start'] + (tss_offset[1] - tss_offset[0])
    # Score column has to be provided for a valid bed file
    df['Score'] = 0
    # Order promoters by chromosome and start
    df.sort_values(['Chromosome', 'Start'], inplace=True)
    
    # Columns to be saved into a file
    columns = ['Chromosome', 'Start', 'End', 'Transcript stable ID', 
        'Score', 'Strand']
    print (f'Saving data to {file_path}...')
    # Save the file
    df[columns].to_csv(file_path, sep='\t', header=False, index=False)
    print ()
    if keep_ensembl:
        # Keep the Ensembl DataFrame in the return dictionary
        answer['ensembl'] = df[['Gene stable ID', 'Transcript stable ID', 
            'Gene name', 'Gene type']]

    return answer


def load_ensembl_from_biomart(
    file_path: FILE_LIKE
) -> Dict[str, Union[str, pd.DataFrame]]:
    """
    Generates the Ensembl file which maps transcripts to genes and 
    stores gene names and types from the data retrieved from the Ensembl
    BioMart server. Returns the database version and the file
    content as a DataFrame.

    Parameters
    ----------
    file_path : FILE_LIKE
        The path to where the resulting file should be saved

    Returns
    -------
    Dict[str, Union[str, pd.DataFrame]]
        A dictionary containing the version of the database used and
        the Ensembl DataFrame
    """
    
    answer = {}

    # Select the right dataset from BioMart
    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    # Attributes to retrieve
    attributes = ['ensembl_transcript_id', 'ensembl_gene_id', 
        'external_gene_name', 'gene_biotype']
    print ('Retrieving response to query...')
    # Submit and retrieve the response
    response = ensembl.search({'attributes': attributes}, header=1)
    buffer = download_with_progress(response)

    # Save the database version into the dictionary
    answer['version'] = ensembl.display_name
    # Convert the response into a DataFrame
    df = pd.read_csv(buffer, sep='\t')

    # Save the file
    df.to_csv(file_path, sep='\t', index=False)
    # Keep the DataFrame in the return dictionary
    answer['ensembl'] = df

    return answer


def download_with_progress(
    url: Union[str, requests.models.Response],
    file_path: Optional[FILE_LIKE] = None,
    desc: str = 'response'
) -> Optional[BytesIO]:
    """
    Downloads from a given URL or retrieves a response to a given 
    request while providing a progress bar.

    Parameters
    ----------
    url : Union[str, requests.models.Response]
        The URL or response to be processed
    file_path : Optional[FILE_LIKE], optional
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
        # Loop until the results are available or max iterations exceeded
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


def get_ensembl_assembly(
) -> str:
    """
    Returns the simple synonym of the genome assembly used by the
    Ensembl server (e.g. hg38).

    Returns
    -------
    str
        A simple synonym of the genome assembly used by Ensembl
    """
    
    # Select the Ensembl dataset from BioMart and get the display name
    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    display_name = ensembl.display_name
    # Isolate the version from the bracket
    version_string = display_name.split('(')[-1].split(')')[0]
    # Remove the update part
    version_major = version_string.split('.')[0]

    # Return the simplified synonym (e.g. hg38 instead of GRCh38):
    return ASSEMBLY_SYNONYM[version_major]


def get_chromosome_mapping(
    assembly: str
) -> Tuple[pd.Series, pd.Series]:
    """
    Returns a tuple with two pandas Series which can be used to map
    Ensembl chromosome names to UCSC (first Series) and vice versa
    (second Series) for a provided genome assembly. If it is not
    recognised, a default mapping valid from the main chromosomes 
    (autosomes + X, Y, MT) is returned.

    Parameters
    ----------
    assembly : str
        The assembly for which to provide the mapping (e.g. hg38)

    Returns
    -------
    Tuple[pd.Series, pd.Series]
        A tuple of two pandas Series, providing chromosome name mapping 
        from Ensembl to UCSC (first one) and vice versa (second one)
    """
    
    if assembly[:2] == 'hg':
        # The mapping can be retrieved from a chromAlias.tsv file
        print (f'Retrieving chromosome name mapping for {assembly}...')
        f = gzip.open(download_with_progress(HG_CHROMOSOME_URL.format(
            genome_assembly=assembly)))
        header_fields = ['alt', 'ucsc', 'notes']
        chrom_df = pd.read_csv(f, sep='\t', names=header_fields)
        # This mapping is unambiguous even if things other than Ensembl
        # are included in the index
        ens_to_ucsc = chrom_df.set_index('alt')['ucsc']
        # For the mapping the other way, we filter out only Ensembl to
        # prevent multiple values for a single index
        chrom_df_filt = chrom_df[chrom_df['notes'].apply(lambda x:
            'ensembl' in x)]
        ucsc_to_ens = chrom_df_filt.set_index('ucsc')['alt']
    else:
        # Default mapping for the 22 autosomal chromosomes + X, Y, MT
        print ('No chromosome name mapping available for the assembly', 
            assembly)
        print ('Using the default mapping')
        ens_to_ucsc = DEFAULT_MAPPING
        # Simple inversion as there is no duplication here
        ucsc_to_ens = pd.Series(ens_to_ucsc.index, index=ens_to_ucsc.values)
    
    return (ens_to_ucsc, ucsc_to_ens)