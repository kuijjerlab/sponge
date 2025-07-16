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
import bioframe
import time

import pandas as pd

from pathlib import Path
from typing import Iterable, Optional

from sponge.modules.utils import iterate_chromosomes, iterate_motifs

### Class definition ###
class MatchFilter:
    # Methods
    def __init__(
        self,
        tfbs_path: Optional[Path],
        regions_path: Path,
        genome_assembly: Optional[str] = None,
        jaspar_release: Optional[str] = None,
        motif_url: Optional[str] = None,
    ):
        """
        Class which filters all the predicted TF binding sites to the
        regions of interest.

        Parameters
        ----------
        tfbs_path : Optional[Path]
            Path to the file containing the TF binding sites or None
            to use on the fly processing instead, in which case
            genome_assembly, jaspar_release and motif_url cannot be None
        regions_path : Path
            Path to the file containing the regions of interest
        genome_assembly : Optional[str], optional
            Genome assembly to use, only used if tfbs_path is None,
            by default None
        jaspar_release : Optional[str], optional
            Version of the JASPAR release to use, only used if tfbs_path
            is None, by default None
        motif_url : Optional[str], optional
            URL to retrieve the TF-specific binding sites from, with
            year and genome_assembly to be formatted in, only used if
            tfbs_path is None, by default None

        Raises
        ------
        ValueError
            If tfbs_path and either genome_assembly, jaspar_release
            or motif_url are None
        """

        if tfbs_path is None:
            error_str = ('At least one of tfbs_path and {parameter} '
                'must be specified')
            if genome_assembly is None:
                raise ValueError(error_str.format(parameter='genome_assembly'))
            if jaspar_release is None:
                raise ValueError(error_str.format(parameter='jaspar_release'))
            if motif_url is None:
                raise ValueError(error_str.format(parameter='motif_url'))

        self.tfbs_path = tfbs_path
        self.regions_path = regions_path
        self.assembly = genome_assembly
        self.jaspar_release = jaspar_release
        self.url = motif_url
        # To be overwritten once retrieved
        self.all_edges = None
        # Overwritten by registering with a VersionLogger instance
        self.version_logger = None


    def filter_matches(
        self,
        chromosomes: Iterable[str],
        matrix_ids: Iterable[str],
        tf_names: Optional[Iterable[str]] = None,
        score_threshold: float = 400,
        n_processes: int = 1,
    ) -> None:
        """
        Filters the binding sites in the bigbed file or on the fly
        downloaded JASPAR tsv files to select only those of the TFs
        of interest in the regions of interest on given chromosomes,
        subject to a score threshold. Stores the result internally.

        Parameters
        ----------
        chromosomes : Iterable[str]
            Which chromosomes to use
        matrix_ids : Iterable[str]
            Which TF matrix IDs to use
        tf_names : Optional[Iterable[str]], optional
            Iterable of TF names matching the matrix IDs, only required
            if on the fly download will be used, by default None
        score_threshold : float, optional
            Minimal score of a match for it to be included in the
            prior, by default 400
        n_processes : int, optional
            Number of processes to run in parallel, by default 1
        """

        print ('\n--- Filtering binding sites in the regions of interest ---')

        print ('Loading the regions file...')
        df_full = bioframe.read_table(self.regions_path, header=0)
        df_full.set_index('Transcript stable ID', inplace=True)

        start_time = time.time()
        if self.tfbs_path is None:
            results_list = iterate_motifs(
                self.url,
                df_full,
                chromosomes,
                matrix_ids,
                tf_names,
                self.assembly,
                self.jaspar_release,
                score_threshold,
                n_processes,
            )
            self.write_retrieved('tfbs_file', self.jaspar_release)
        else:
            results_list = iterate_chromosomes(
                self.tfbs_path,
                df_full,
                chromosomes,
                matrix_ids,
                score_threshold,
                n_processes,
            )

        elapsed = time.time() - start_time

        print (f'\nTotal time: {elapsed // 60:n} m {elapsed % 60:.2f} s')

        # Save the final results, ignoring the index makes this fast
        # The index is irrelevant
        self.all_edges = pd.concat(results_list, ignore_index=True)


    def get_all_edges(
        self,
    ) -> pd.DataFrame:
        """
        Get the pandas DataFrame of all the filtered edges, can be None
        if it wasn't created yet.

        Returns
        -------
        pd.DataFrame
            DataFrame containing all the filtered edges
        """

        return self.all_edges

    # Placeholder functions to be replaced with VersionLogger if any
    def write_retrieved(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass