### Imports ###
import os
import yaml

from collections import defaultdict
from pathlib import Path
from typing import Optional
from yaml.error import MarkedYAMLError

### Class definition ###
class VersionLogger:
    # Variables
    _log_filename = 'fingerprint.yaml'

    # Functions
    def __init__(
        self,
        temp_folder: Path,
        log_filename: Optional[str] = None,
    ):

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


    def __del__(
        self,
    ):

        if len(self.data) > 0:
            yaml.safe_dump(dict(self.data), open(self.log_file, 'w'))


    def __getitem__(
        self,
        key: str,
    ) -> dict:

        return self.data[key]


    def __setitem__(
        self,
        key: str,
        val: dict,
    ) -> None:

        self.data[key] = val


    def __delitem__(
        self,
        key: str,
    ) -> None:

        del self.data[key]

    
    def write_provided(
        self,
        key: str,
    ) -> None:
        
        pass


    def update_cached(
        self,
        key: str,
    ) -> None:
        
        pass


    def write_retrieved(
        self,
        key: str,
        version: str,
    ) -> None:
        
        pass