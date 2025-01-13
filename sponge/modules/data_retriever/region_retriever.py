### Imports ###
import os

import pandas as pd

from collections import defaultdict
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.data_retriever.file_retriever import FileRetriever
from sponge.modules.utils import get_ensembl_version, retrieve_ensembl_data, \
    get_chromosome_mapping
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class RegionRetriever(FileRetriever):
    _default_filename = 'regions.tsv'

    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigReader,
        user_config: ConfigReader,
        version_logger: VersionLogger,
    ):

        path_to_file = None
        if user_config.exists(['region', 'region_file']):
            path_to_file = user_config.get_value(['region', 'region_file'])
        temp_filename = os.path.join(temp_folder, self._default_filename)

        self.xml = core_config['url']['region']['xml']
        self.rest = core_config['url']['region']['rest']
        self.mapping = core_config['default_mapping']
        self.assembly = user_config['genome_assembly']
        self.mapping_url = core_config['url']['chrom_mapping']

        self.settings = user_config['region']
        if 'chromosomes' not in self.settings:
            self.settings['chromosomes'] = core_config['default_chromosomes']

        super().__init__(
            key='region_file',
            temp_filename=temp_filename,
            version_logger=version_logger,
            path_to_file=path_to_file,
        )


    def retrieve_file(
        self,
        temp_filename: Path,
    ) -> str:

        # Attributes to retrieve
        attributes = ['ensembl_transcript_id', 'transcript_gencode_basic',
            'chromosome_name', 'transcription_start_site', 'strand',
            'ensembl_gene_id', 'external_gene_name', 'gene_biotype']
        # Submit and retrieve the response
        buffer = retrieve_ensembl_data('hsapiens_gene_ensembl', attributes,
            self.xml)

        # Dictionary of types for conversion from the response, default strings
        dtype_dict = defaultdict(lambda: str)
        # Change the types that are not strings but integers
        dtype_dict['Transcription start site (TSS)'] = int
        dtype_dict['Strand'] = int
        # Convert the response into a DataFrame
        df = pd.read_csv(buffer, sep='\t', dtype=dtype_dict)

        print ('Filtering and modifying dataframe...')
        if self.settings['filter_basic']:
            # Filter only for GENCODE basic
            df = df[df['GENCODE basic annotation'] == 'GENCODE basic'].copy()
            df.drop(columns='GENCODE basic annotation', inplace=True)
        # Convert chromosome names to match with other inputs
        # Attempt to retrieve mapping
        mapping = get_chromosome_mapping(self.assembly, self.mapping_url)
        if mapping is not None:
            # Managed to retrieve it, overwrite the default
            self.mapping = mapping
        df['Chromosome'] = df['Chromosome/scaffold name'].apply(lambda x:
            self.mapping[x] if x in self.mapping else None)
        not_mapped = df['Chromosome'].isna().sum()
        if not_mapped > 0:
            print (f'Discarding {not_mapped} regions on unmapped chromosomes.')
            df.dropna(subset=['Chromosome'], inplace=True)
        chromosomes = self.settings['chromosomes']
        if chromosomes is not None:
            # Filter only for selected chromosomes
            df = df[df['Chromosome'].isin(chromosomes)]
        # Convert strand to +/-
        df['Strand'] = df['Strand'].apply(lambda x: '+' if x > 0 else '-')
        # Calculate the start based on the given offset from TSS
        # The calculation is strand dependent
        tss_offset = self.settings['tss_offset']
        df['Start'] = df.apply(lambda row:
            row['Transcription start site (TSS)'] + tss_offset[0]
            if row['Strand'] == '+'
            else row['Transcription start site (TSS)'] - tss_offset[1],
            axis=1)
        # End is always greater than start, this way it is strand independent
        df['End'] = df['Start'] + (tss_offset[1] - tss_offset[0])
        # Order promoters by chromosome and start
        df.sort_values(['Chromosome', 'Start'], inplace=True)

        # Columns to be saved into a file
        columns = ['Chromosome', 'Start', 'End', 'Transcript stable ID',
            'Gene stable ID', 'Gene name', 'Gene type']
        print (f'Saving data to: {temp_filename}')
        # Save the file
        self.df = df[columns]
        self.df.to_csv(temp_filename, sep='\t', index=False)
        print ()

        return get_ensembl_version(self.rest)