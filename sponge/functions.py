### Imports ###
import pandas as pd
import bioframe
import requests
import os
import time

from Bio.motifs.jaspar import Motif

from typing import Optional, Union, Iterable, Tuple

from biomart import BiomartServer

from io import BytesIO

from collections import defaultdict

from math import log2

from tqdm import tqdm

from io import BytesIO

### Data ###
FILE_DF = pd.DataFrame(
    {'description': ['homologene', 'promoter', 'jaspar_bigbed', 'ensemble'],
     'name': ['homologene.tsv', 'promoters.bed', 'JASPAR.bb', 'ensemble.tsv'],
     'url': ['https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data',
             None,
             'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/{year}'
             '/JASPAR{year}_hg38.bb',
             None],
     'eval': [None,
              'self.load_promoters_from_biomart(**options)',
              None,
              'self.load_ensemble_from_biomart(**options)']}
).set_index('description')

ENSEMBL_URL = 'http://www.ensembl.org/biomart'
MAPPING_URL = 'https://rest.uniprot.org/idmapping/'

### Functions ###
def download_with_progress(
    url: Union[str, requests.models.Response],
    file_path: Optional[str] = None,
    desc: str = 'response'
) -> Optional[BytesIO]:
    
    if type(url) == str:
        request = requests.get(url, stream=True)
    else:
        request = url
    total = int(request.headers.get('content-length', 0))          
    if file_path is None:
        stream = BytesIO()
    else:
        stream = open(file_path, 'wb')
        desc = file_path
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
    
    data = {'ids': ids, 'from': from_db, 'to': to_db}
    data.update(kwargs)
    uniprot_request = requests.post(MAPPING_URL + 'run', data)
    uniprot_reply = uniprot_request.json()
    if 'jobId' in uniprot_reply:
        job_id = uniprot_reply['jobId']
    else:
        print ('Unable to retrieve a job ID from UniProt')
        print ('The following reply was received instead:')
        for message in uniprot_reply['messages']:
            print (message)
        raise RuntimeError()
    while True:
        uniprot_status = requests.get(MAPPING_URL + f'status/{job_id}')
        if 'results' in uniprot_status.json():
            break
        time.sleep(0.5)
    uniprot_results = requests.get(MAPPING_URL + f'stream/{job_id}')
    results_df = pd.DataFrame(uniprot_results.json()['results'])
    results_df.drop_duplicates(subset='from', inplace=True)
    return results_df


def process_ensemble_df(
    df: pd.DataFrame
) -> None:

    pass


def load_promoters_from_biomart(
    file_path: str,
    filter_basic: bool = True,
    chromosomes: Iterable[str] = [str(i) for i in range(1,23)] + 
        ['MT', 'X', 'Y'],
    save_ensemble: bool = True
) -> None:

    bm_server = BiomartServer(ENSEMBL_URL)
    ensembl = bm_server.datasets['hsapiens_gene_ensembl']
    attributes = ['ensembl_transcript_id', 'transcript_gencode_basic', 
        'chromosome_name', 'transcription_start_site', 'strand']
    if save_ensemble:
        attributes += ['ensembl_gene_id', 'external_gene_name', 'gene_biotype']
    print ('Retrieving response to query...')
    dtype_dict = defaultdict(lambda: str)
    dtype_dict['Transcription start site (TSS)'] = int
    dtype_dict['Strand'] = int
    response = ensembl.search({'attributes': attributes}, header=1)
    buffer = download_with_progress(response)
    df = pd.read_csv(buffer, sep='\t', dtype=dtype_dict)
    print ('Filtering and modifying dataframe...')
    if filter_basic:
        df = df[df['GENCODE basic annotation'] == 'GENCODE basic'].copy()
        df.drop(columns='GENCODE basic annotation', inplace=True)
    if chromosomes is not None:
        df = df[df['Chromosome/scaffold name'].isin(chromosomes)]
    df['Chromosome'] = df['Chromosome/scaffold name'].apply(lambda x: 'chrM' 
        if x == 'MT' else f'chr{x}')
    df['Strand'] = df['Strand'].apply(lambda x: '+' if x > 0 else '-')
    df['Start'] = df.apply(lambda row: 
        row['Transcription start site (TSS)'] - 750 if row['Strand'] == '+' 
        else row['Transcription start site (TSS)'] - 250, axis=1)
    df['End'] = df['Start'] + 1000
    df['Score'] = 0
    df.sort_values(['Chromosome', 'Start'], inplace=True)
    columns = ['Chromosome', 'Start', 'End', 'Transcript stable ID', 'Score', 
        'Strand']
    print (f'Saving data to {file_path}...')
    df[columns].to_csv(file_path, sep='\t', header=False, index=False)
    print ()
    if save_ensemble:
        columns = ['Gene stable ID', 'Transcript stable ID', 'Gene name', 
            'Gene type']
        process_ensemble_df(df[columns])


def load_ensemble_from_biomart(
    file_path: str
) -> None:
    
    pass


def retrieve_file(
    description: str,
    temp_folder: str,
    prompt: bool = True,
    jaspar_release: Optional[str] = None
) -> str:

    
    if description not in FILE_DF.index:
        print (f'File description not recognised: {description}')
        return
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    file_name = FILE_DF.loc[description, 'name']
    file_path = os.path.join(temp_folder, file_name)
    if os.path.exists(file_path):
        print (f'Using cached file {file_path}')
        print ()
    else:
        print (f'File {file_name} not found in directory {temp_folder}')
        if prompt:
            key = None
            positive = ['y', 'yes', 'hell yeah']
            negative = ['n', 'no', 'nope']
            while key is None or key.lower() not in positive + negative:
                if key is not None:
                    print (f'Input not recognised: {key}')
                print ('Do you want to download it? Y/N', flush=True)
                key = input()
                print (key)
            print ()
            if key.lower() in negative:
                return None
        to_request = FILE_DF.loc[description, 'url']
        if to_request is None:
            # These options are not unused: they are passed to the
            # evaluated function call as kwargs
            options = {'file_path': file_path}
            eval(FILE_DF.loc[description, 'eval'])
        else:
            if description == 'jaspar_bigbed':
                if jaspar_release is None:
                    raise ValueError('The release of jaspar has to be '
                        'specified in order to retrieve the bigbed file')
                to_request = to_request.format(year=jaspar_release[-4:])
            print (f'Downloading data into {file_path}...')
            download_with_progress(to_request, file_path)    

    return file_path


def filter_edges(
    bb_ref: str, 
    bed_df: pd.DataFrame, 
    motif_list: Iterable[str], 
    chrom: str, 
    start_ind: int, 
    final_ind: int, 
    score_threshold: float = 400
) -> pd.DataFrame:
    
    df = pd.DataFrame()
    for transcript in bed_df.index[start_ind:final_ind]:
        start,end = bed_df.loc[transcript][['start', 'end']]
        motifs = bioframe.read_bigbed(bb_ref, chrom, start=start, end=end)
        max_scores = motifs[motifs['TFName'].isin(motif_list)].groupby(
            'TFName')[['score', 'name']].max()
        max_scores = max_scores[max_scores['score'] >= score_threshold]
        max_scores.reset_index(inplace=True)
        max_scores['transcript'] = transcript
        df = pd.concat((df, max_scores), ignore_index=True, copy=False)
    return df


def filter_edges_helper(
    input_tuple: Tuple[str, pd.DataFrame, Iterable[str], str, int, int]
) -> pd.DataFrame:

    return filter_edges(*input_tuple)


def plogp(
    x: float
) -> float:   
    """
    Returns x*log2(x) for a number, handles the 0 case properly.
    """

    if x == 0:
        return 0
    else:
        return x*log2(x)


def calculate_ic(
    motif: Motif
) -> float:
    """
    Calculates the information content for a given motif, assumes equal 
    ACGT distribution.
    """

    df = pd.DataFrame(motif.pwm)
    df['IC'] = df.apply(lambda x: 2 + sum([plogp(x[y]) for y in 
        ['A', 'C', 'G', 'T']]), axis=1)
    
    return df['IC'].sum()


def adjust_gene_name(
    gene: str
) -> str:
    
    return gene[:-2] + gene[-2:].lower()