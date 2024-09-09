import pandas as pd

# The links for downloading files or alternatively the names of functions that
# should be executed in order to retrieve them
FILE_DF = pd.DataFrame(
    {'description': ['promoter', 'jaspar_bigbed', 'ensembl'],
     'name': ['promoters.bed', 'JASPAR.bb', 'ensembl.tsv'],
     'url': [None,
             ['https://frigg.uio.no/JASPAR/JASPAR_TFBSs/{year}'
             '/JASPAR{year}_{genome_assembly}.bb',
             'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/{year}'
             '/JASPAR{year}_{genome_assembly}.bb'],
             None],
     'eval': ['load_promoters_from_biomart(**options)',
              None,
              'load_ensembl_from_biomart(**options)']}
).set_index('description')

# URLs to websites for downloads, should only be provided here and referenced
# by their variable names for ease of future updates
MOTIF_URL = ['https://frigg.uio.no/JASPAR/JASPAR_TFBSs/{year}/'
    '{genome_assembly}/',
    'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/{year}/'
    '{genome_assembly}/']
ENSEMBL_URL = 'http://www.ensembl.org/biomart'
MAPPING_URL = 'https://rest.uniprot.org/idmapping/'
STRING_URL = 'https://string-db.org/api/tsv/'
NCBI_URL = 'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/'
HG_CHROMOSOME_URL = ('https://hgdownload.soe.ucsc.edu/goldenPath/'
    '{genome_assembly}/database/chromAlias.txt.gz')

# Synonyms for genome assembly versions
ASSEMBLY_SYNONYM = {'GRCh38': 'hg38', 'GRCh37': 'hg19', 'T2T-CHM13v2.0': 'hs1'}

# Default chromosomes to use
DEFAULT_CHROMOSOMES = [f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']]

# Default chromosome name mapping from Ensembl to UCSC
index = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
values = ['chrM' if i == 'MT' else f'chr{i}' for i in index]
DEFAULT_MAPPING = pd.Series(values, index=index)