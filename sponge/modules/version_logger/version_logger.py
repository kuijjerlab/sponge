### Imports ###
import os
import yaml

from pathlib import Path

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

        if os.path.exists(self.log_file):
            self.data = yaml.safe_load(open(self.log_file))
        else:
            self.data = {}


    def __del__(
        self
    ):

        yaml.safe_dump(self.data, open(self.log_file, 'w'))