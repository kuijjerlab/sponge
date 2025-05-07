### Imports ###
import os

import pandas as pd

from collections import defaultdict
from pathlib import Path
from typing import List, Mapping

from sponge.modules.data_retriever.file_retriever import FileRetriever
from sponge.modules.utils import get_ensembl_version, retrieve_ensembl_data, \
    get_chromosome_mapping

### Class definition ###
class RegionRetriever(FileRetriever):
    # Variables
    _default_filename = 'regions.tsv'

    # Methods
    def __init__(
        self,
        temp_folder: Path,
        xml_url: str,
        rest_url: str,
        mapping_url: str,
        region_settings: dict,
        genome_assembly: str,
        default_mapping: Mapping,
        default_chromosomes: List[str],
    ):
        """
        Class which retrieves suitable regions of interest given the
        provided settings.

        Parameters
        ----------
        temp_folder : Path
            Folder to save the retrieved file to
        xml_url : str
            URL to retrieve Ensembl xml file from (likely BioMart)
        rest_url : str
            URL for REST interface to Ensembl
        mapping_url : str
            URL to retrieve the mapping from, with genome_assembly to
            be formatted in
        region_settings : dict
            Region relevant settings from the user configuration,
            specifically containing region_file (the path to the
            provided region file), filter_basic, use_all_chromosomes,
            chromosomes and tss_offset
        genome_assembly : str
            Genome assembly to use
        default_mapping : Mapping
            Mapping between Ensembl and UCSC chromosome names which
            will be used by default
        default_chromosomes : List[str]
            List of chromosomes which will be used by default
        """

        temp_filename = os.path.join(temp_folder, self._default_filename)

        self.xml = xml_url
        self.rest = rest_url
        self.mapping_url = mapping_url
        self.assembly = genome_assembly
        self.mapping = default_mapping

        if region_settings['chromosomes'] is None:
            region_settings['chromosomes'] = default_chromosomes
        self.settings = region_settings

        super().__init__(
            key='region_file',
            temp_filename=temp_filename,
            path_to_file=region_settings['region_file'],
        )


    def _retrieve_region(
        self,
    ) -> str:
        """
        Retrieves the region file from Ensembl.

        Returns
        -------
        str
            Version of the Ensembl database corresponding to the file
        """

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
        # Filter for the desired chromosomes
        if self.settings['use_all_chromosomes']:
            # Not used here but SPONGE may use it downstream
            self.settings['chromosomes'] = list(df['Chromosome'].unique())
        else:
            df = df[df['Chromosome'].isin(self.settings['chromosomes'])]
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
        print (f'Saving data to: {self.temp_filename}')
        # Save the file
        self.df = df[columns]
        self.df.to_csv(self.temp_filename, sep='\t', index=False)

        return get_ensembl_version(self.rest)


    def retrieve_file(
        self,
    ) -> None:
        """
        Retrieves the region file.
        """

        super().retrieve_file(self._retrieve_region)