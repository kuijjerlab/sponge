import pandas as pd

# The links for downloading files or alternatively the names of functions that 
# should be executed in order to retrieve them
FILE_DF = pd.DataFrame(
    {'description': ['homologene', 'promoter', 'jaspar_bigbed', 'ensembl'],
     'name': ['homologene.tsv', 'promoters.bed', 'JASPAR.bb', 'ensembl.tsv'],
     'url': ['https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data',
             None,
             'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/{year}'
             '/JASPAR{year}_{genome_assembly}.bb',
             None],
     'eval': [None,
              'load_promoters_from_biomart(**options)',
              None,
              'load_ensembl_from_biomart(**options)']}
).set_index('description')

# URLs to websites for downloads, should only be provided here and referenced
# by their variable names for ease of future updates
ENSEMBL_URL = 'http://www.ensembl.org/biomart'
MAPPING_URL = 'https://rest.uniprot.org/idmapping/'
STRING_URL = 'https://string-db.org/api/tsv/'
HG_CHROMOSOME_URL = ('https://hgdownload.soe.ucsc.edu/goldenPath/'
    '{genome_assembly}/database/chromAlias.txt.gz')

# Synonyms for genome assembly versions
ASSEMBLY_SYNONYM = {'GRCh38': 'hg38', 'GRCh37': 'hg19', 'T2T-CHM13v2.0': 'hs1'}

# Default chromosome name mapping from Ensembl to UCSC
index = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']
values = ['chrM' if i == 'MT' else f'chr{i}' for i in index]
DEFAULT_MAPPING = pd.Series(values, index=index)

# TODO: unipressed module to replace Uniprot calls?