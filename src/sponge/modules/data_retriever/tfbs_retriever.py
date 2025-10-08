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
import os

from pathlib import Path
from typing import List, Union

from sponge.modules.data_retriever.file_retriever import FileRetriever
from sponge.modules.utils import download_with_progress

### Class definition ###
class TFBSRetriever(FileRetriever):
    # Variables
    _default_filename = 'tfbs.bb'

    # Methods
    def __init__(
        self,
        temp_folder: Path,
        motif_url: Union[str, List[str]],
        motif_settings: dict,
        genome_assembly: dict,
        on_the_fly_processing: bool = True,
    ):
        """
        Class which retrieves the TFBS file from JASPAR given the
        provided settings.

        Parameters
        ----------
        temp_folder : Path
            Folder to save the retrieved file to
        motif_url : Union[str, List[str]]
            URL or list of URLS to try and retrieve the TFBS file from,
            with year and genome_assembly to be formatted in
        motif_settings : dict
            Motif relevant settings from the user configuration,
            specifically containing jaspar_release and tfbs_file (the
            path to the provided TFBS file)
        genome_assembly : dict
            Genome assembly to use
        on_the_fly_processing : bool, optional
            Whether to process the TFBS with on the fly download of
            individual TF files which would mean the TFBS file will not
            be downloaded, by default True
        """

        self.on_the_fly = on_the_fly_processing
        self.jaspar_release = motif_settings['jaspar_release']
        self.genome_assembly = genome_assembly
        self.urls = motif_url

        temp_filename = os.path.join(temp_folder, self._default_filename)

        super().__init__(
            key='tfbs_file',
            temp_filename=temp_filename,
            path_to_file=motif_settings['tfbs_file'],
        )


    def _retrieve_tfbs(
        self,
    ) -> str:
        """
        Retrieves the TFBS file by downloading it.

        Returns
        -------
        str
            Version of the JASPAR release corresponding to the file
        """

        year = self.jaspar_release[-4:]
        urls_to_try = [url.format(year=year,
            genome_assembly=self.genome_assembly) for url in self.urls]
        download_with_progress(urls_to_try, self.temp_filename)

        return self.jaspar_release


    def retrieve_file(
        self,
    ) -> None:
        """
        Retrieves the TFBS file or skips that if it is not required.
        """

        if self.on_the_fly:
            print ('Retrieval of tfbs_file is skipped as on the fly '
                'processing was requested.')
        else:
            super().retrieve_file(self._retrieve_tfbs)