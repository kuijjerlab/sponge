### Imports ###
import os
import yaml

from collections import defaultdict
from pathlib import Path
from yaml.error import MarkedYAMLError

### Class definition ###
class VersionLogger:
    # Variables
    _log_filename = 'fingerprint.yaml'

    # Functions
    def __init__(
        self,
        temp_folder: Path,
    ):

        os.makedirs(temp_folder, exist_ok=True)
        self.log_file = os.path.join(temp_folder, self._log_filename)

        self.data = defaultdict(dict)
        if os.path.exists(self.log_file):
            try:
                for k,v in yaml.safe_load(open(self.log_file)).items():
                    self.data[k] = v
            except AttributeError:
                # Most likely means an empty log file, ignore
                pass
            except MarkedYAMLError:
                print ('There seems to be an issue with the fingeprint file. '
                    f'We recommend deleting the temporary folder {temp_folder}'
                    ' to fix the issue.')


    def __del__(
        self
    ):

        temp_dict = { k: v for k,v in self.data.items() }
        yaml.safe_dump(temp_dict, open(self.log_file, 'w'))


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