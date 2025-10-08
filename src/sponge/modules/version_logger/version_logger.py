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
import datetime
import os
import time
import yaml

from collections import defaultdict
from pathlib import Path
from typing import Any, Optional
from yaml.error import MarkedYAMLError

### Class definition ###
class VersionLogger:
    # Variables
    _log_filename = 'fingerprint.yaml'

    # Methods
    def __init__(
        self,
        temp_folder: Path,
        log_filename: Optional[str] = None,
    ):
        """
        Logs the accession times and versions of databases and files
        and maintains a file with the records. 

        Parameters
        ----------
        temp_folder : Path
            Folder where the fingerprint file should be maintained
        log_filename : Optional[str], optional
            Name of the fingerprint file or None to use the default,
            by default None
        """

        os.makedirs(temp_folder, exist_ok=True)
        if log_filename is None:
            log_filename = self._log_filename
        self.log_file = os.path.join(temp_folder, log_filename)

        self.data = defaultdict(dict)
        if os.path.exists(self.log_file):
            try:
                self.data = defaultdict(dict,
                    yaml.safe_load(open(self.log_file)))
            except TypeError:
                # Most likely means an empty log file, ignore
                pass
            except MarkedYAMLError:
                print(
                    'There seems to be an issue with the fingeprint file. '
                    f'We recommend deleting the temporary folder {temp_folder}'
                    ' to fix the issue.'
                )

    # Overwritten built-in methods
    def __del__(
        self,
    ):
        """
        Saves the records into a file when the object is deleted,
        assuming there is any.
        """

        if len(self.data) > 0:
            yaml.safe_dump(dict(self.data),
                open(self.log_file, 'w', encoding='utf-8'))


    def __getitem__(
        self,
        key: str,
    ) -> dict:
        """
        Retrieves the values for a single record.

        Parameters
        ----------
        key : str
            Key describing the record

        Returns
        -------
        dict
            Values associated with the record
        """

        return self.data[key]


    def __setitem__(
        self,
        key: str,
        val: dict,
    ) -> None:
        """
        Sets the value of a single record.

        Parameters
        ----------
        key : str
            Key describing the record
        val : dict
            Values to be stored in it
        """

        self.data[key] = val


    def __delitem__(
        self,
        key: str,
    ) -> None:
        """
        Deletes a single record.

        Parameters
        ----------
        key : str
            Key describing the record
        """

        del self.data[key]


    def __contains__(
        self,
        key: str,
    ) -> bool:
        """
        Checks whether a specified record exists.

        Parameters
        ----------
        key : str
            Key describing the record

        Returns
        -------
        bool
            Whether the record exists
        """

        return key in self.data

    # Rest of the methods
    def _reset_entry(
        self,
        key: str,
    ) -> None:
        """
        Resets the entry for the given record to the default, which
        means it deletes the already present data and sets retrieval
        time to current time.

        Parameters
        ----------
        key : str
            Key describing the record
        """

        # Remove if present
        if key in self.data:
            del self.data[key]

        raw_dt = datetime.datetime.fromtimestamp(time.time())
        self.data[key]['datetime'] = raw_dt.replace(microsecond=0)


    def write_provided(
        self,
        key: str,
    ) -> None:
        """
        Specifies the given record as provided by the user at the
        current time.

        Parameters
        ----------
        key : str
            Key describing the record
        """

        self._reset_entry(key)

        self.data[key]['provided'] = True
        self.data[key]['version'] = 'unknown'


    def write_default(
        self,
        key: str,
    ) -> None:
        """
        Specifies the given record as default, which means neither
        the time it was provided nor the version are known.

        Parameters
        ----------
        key : str
            Key describing the record
        """

        self._reset_entry(key)
        # Datetime is set here again
        self.data[key]['datetime'] = 'unknown'
        self.data[key]['version'] = 'unknown'


    def update_cached(
        self,
        key: str,
    ) -> None:
        """
        Updates the given record to specify it was retrieved from cache.

        Parameters
        ----------
        key : str
            Key describing the record
        """

        # Only update cache label, don't change anything else
        self.data[key]['cached'] = True


    def write_retrieved(
        self,
        key: str,
        version: str,
    ) -> None:
        """
        Specifies the given record as retrieved at the current time with
        the given version.

        Parameters
        ----------
        key : str
            Key describing the record
        version : str
            Version of the record
        """

        self._reset_entry(key)

        self.data[key]['version'] = version


    def register_class(
        self,
        target_class: type,
    ) -> None:
        """
        Registers a class instance with this instance, giving it access
        to the logging functionality.

        Parameters
        ----------
        target_class : type
            Class instance to be modified
        """

        # Include a reference to this instance in the registered class
        target_class.version_logger = self
        # List of functions to be overwritten in the registered class
        function_names = ['write_provided', 'write_default', 'update_cached',
            'write_retrieved']
        for fn in function_names:
            # Replace with a call to the internal version logger
           
            setattr(target_class, fn, getattr(self, fn))