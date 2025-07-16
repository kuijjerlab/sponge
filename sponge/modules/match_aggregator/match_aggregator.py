# Copyright (C) 2025 Ladislav Hovan <ladislav.hovan@ncmbm.uio.no>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public License as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public License along
# with this library. If not, see <https://www.gnu.org/licenses/>.

### Imports ###
import numpy as np
import pandas as pd

from pathlib import Path
from typing import Mapping

### Class definition ###
class MatchAggregator:
    # Methods
    def __init__(
        self,
        edges: pd.DataFrame,
        regions_file: Path,
        homolog_mapping: Mapping[str, str],
    ):
        """
        Class that aggregates transcript matches to genes and can also
        filter the genes and convert IDs to names.

        Parameters
        ----------
        edges : pd.DataFrame
            DataFrame containing the edges of the motif prior
        regions_file : Path
            Path to the file describing the regions of interest
        homolog_mapping : Mapping[str, str]
            Mapping of the TF homologs in other species to the
            species of interest
        """

        self.initial_edges = edges
        self.regions = pd.read_csv(regions_file, sep='\t')
        self.homolog_mapping = homolog_mapping
        # To be overwritten once retrieved
        self.edges = None


    def aggregate_matches(
        self,
        use_gene_names: bool = True,
        protein_coding_only: bool = False,
    ) -> pd.DataFrame:
        """
        Aggregates all the matches corresponding to individual
        transcripts into genes, creating a transcription factor - gene
        matrix. Stores the result internally.

        Parameters
        ----------
        use_gene_names : bool, optional
            Whether to use gene names instead of Ensembl IDs,
            by default True
        protein_coding_only : bool, optional
            Whether to restrict the selection to only protein coding
            genes, by default False
        """

        print ('\n--- Aggregating transcript matches into genes ---')

        # Add the gene name data to the edges previously found
        motif_df = self.initial_edges.join(other=self.regions.set_index(
            'Transcript stable ID'), on='transcript')
        print ('Number of TF - transcript edges:', len(motif_df))
        if protein_coding_only:
            motif_df = motif_df[motif_df['Gene type'] ==
                'protein_coding'].copy()
        # Drop columns that are not required anymore
        motif_df.drop(columns=['Gene type', 'name'], inplace=True)
        # Humanise the TF names
        motif_df['TFName'] = motif_df['TFName'].apply(lambda x:
            self.homolog_mapping[x] if x in self.homolog_mapping else x)
        # Ignore genes without identifiers
        motif_df.dropna(subset=['Gene stable ID'], inplace=True)
        motif_df.sort_values('score', ascending=False, inplace=True)
        # Sometimes edges are identified from multiple transcripts
        motif_df.drop_duplicates(subset=['TFName', 'Gene stable ID'],
            inplace=True)
        print ('Number of TF - gene edges:', len(motif_df))
        if use_gene_names:
            # Names are not unique - filtering needed
            # Fill empty gene names with IDs
            motif_df['Gene name'] = motif_df.apply(lambda x: x['Gene name'] if
                type(x['Gene name']) == str else x['Gene stable ID'], axis=1)
            # Count the number of edges for every name/ID pair
            name_id_matching = motif_df.groupby(
                ['Gene name', 'Gene stable ID'])['Gene name'].count()
            # Use the name for the ID that has the most edges
            id_to_name = {i[1]: i[0] for i in name_id_matching.groupby(
                level=0).idxmax().values}
            # Convert selected gene IDs to names
            motif_df['Gene name'] = motif_df['Gene stable ID'].apply(
                lambda x: id_to_name[x] if x in id_to_name else np.nan)
            # Drop the rest
            motif_df.dropna(subset='Gene name', inplace=True)
            print ('Number of TF - gene edges after name conversion:',
                len(motif_df))

        self.edges = motif_df


    def get_edges(
        self,
    ) -> pd.DataFrame:
        """
        Get the final edges of the prior gene regulatory network,
        can be None if it wasn't generated yet.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the gene regulatory network
        """

        return self.edges