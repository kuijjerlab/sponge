### Imports ###
import datetime
import os
import time
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
            yaml.safe_dump(dict(self.data),
                open(self.log_file, 'w', encoding='utf-8'))


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


    def __contains__(
        self,
        key: str,
    ) -> bool:

        return key in self.data


    def _reset_entry(
        self,
        key: str,
    ) -> None:

        # Remove if present
        if key in self.data:
            del self.data[key]

        raw_dt = datetime.datetime.fromtimestamp(time.time())
        self.data[key]['datetime'] = raw_dt.replace(microsecond=0)


    def write_provided(
        self,
        key: str,
    ) -> None:

        self._reset_entry(key)

        self.data[key]['provided'] = True
        self.data[key]['version'] = 'unknown'


    def write_default(
        self,
        key: str,
    ) -> None:

        self._reset_entry(key)
        # Datetime is set here again
        self.data[key]['datetime'] = 'unknown'
        self.data[key]['version'] = 'unknown'


    def update_cached(
        self,
        key: str,
    ) -> None:

        # Only update cache label, don't change anything else
        self.data[key]['cached'] = True


    def write_retrieved(
        self,
        key: str,
        version: str,
    ) -> None:

        self._reset_entry(key)

        self.data[key]['version'] = version