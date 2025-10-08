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
from typing import Callable, Optional

### Class definition ###
class FileRetriever:
    # Methods
    def __init__(
        self,
        key: str,
        temp_filename: str,
        path_to_file: Optional[Path] = None,
    ):
        """
        Base class for the retrieval of files required by SPONGE.
        Takes care of checking the cache and updating the values
        if registered by a VersionLogger.

        Parameters
        ----------
        key : str
            Descriptor of the file to be retrieved
        temp_filename : str
            Path to where a cached file would be located
        path_to_file : Optional[Path], optional
            Path to the provided file or None if none was provided,
            by default None
        """

        self.key = key
        self.temp_filename = temp_filename
        self.path_to_file = path_to_file
        # Meant to be overwritten after retrieval
        self.actual_path = None
        # Overwritten by registering with a VersionLogger instance
        self.version_logger = None


    def retrieve_file(
        self,
        retrieve_function: Callable,
    ) -> None:
        """
        Attempts to retrieve the specified file from the cache or
        the user-provided path, retrieves it if not found, logs the
        results if registered by a VersionLogger.

        Parameters
        ----------
        retrieve_function : Callable
            Function to be called to retrieve the file if it is not
            found in the cache or provided directly, it should take
            no arguments and return the version

        Raises
        ------
        FileNotFoundError
            If the user-provided path does not exist
        """

        print (f'\n--- Attempting to locate {self.key} ---')

        # Check existence of a user-provided file
        if self.path_to_file is not None:
            print (f'Using a user-provided file: {self.path_to_file}')
            if not os.path.exists(self.path_to_file):
                raise FileNotFoundError('Could not locate file: '
                    f'{self.path_to_file}')
            self.write_provided(self.key)
            self.actual_path = self.path_to_file
        # Check for a cached file
        elif os.path.exists(self.temp_filename):
            if (self.version_logger is not None and
                self.key not in self.version_logger):
                print ('A cached file is present but it is not being tracked '
                    'by the version logger.')
                self.write_default(self.key)
            print ('Reusing a cached file.')
            self.update_cached(self.key)
            self.actual_path = self.temp_filename
        # Retrieve a file
        else:
            print ('Retrieving the file...')
            version = retrieve_function()
            self.write_retrieved(self.key, version)
            self.actual_path = self.temp_filename

    # Placeholder functions to be replaced with VersionLogger if any
    def write_provided(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass


    def write_default(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass


    def update_cached(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass


    def write_retrieved(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass