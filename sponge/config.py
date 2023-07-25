import pandas as pd

# The links for downloading files or alternatively the names of functions that 
# should be executed in order to retrieve them
FILE_DF = pd.DataFrame(
    {'description': ['homologene', 'promoter', 'jaspar_bigbed', 'ensembl'],
     'name': ['homologene.tsv', 'promoters.bed', 'JASPAR.bb', 'ensembl.tsv'],
     'url': ['https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data',
             None,
             'http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/{year}'
             '/JASPAR{year}_hg38.bb',
             None],
     'eval': [None,
              'self.load_promoters_from_biomart(**options)',
              None,
              'self.load_ensembl_from_biomart(**options)']}
).set_index('description')

# URLs to websites for downloads, should only be provided here and referenced
# by their variable names for ease of future updates
ENSEMBL_URL = 'http://www.ensembl.org/biomart'
MAPPING_URL = 'https://rest.uniprot.org/idmapping/'
STRING_URL = 'https://string-db.org/api/tsv/'

# TODO: unipressed module to replace Uniprot calls?