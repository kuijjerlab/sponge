### Imports ###
import os

from pathlib import Path
from typing import Optional

from sponge.modules.version_logger import VersionLogger

### Class definition ###
class FileRetriever:
    def __init__(
        self,
        key: str,
        temp_filename: str,
        version_logger: VersionLogger,
        path_to_file: Optional[Path] = None,
    ):

        print ()
        print (f'--- Attempting to locate {key} file ---')

        # Check existence of a user-provided file
        if path_to_file is not None:
            print (f'Using a user-provided file: {path_to_file}')
            if not os.path.exists(path_to_file):
                raise FileNotFoundError('could not locate file: '
                    f'{path_to_file}')
            version_logger.write_provided(key)
        # Check for a cached file
        elif os.path.exists(temp_filename):
            if key not in version_logger:
                print ('A cached file is present but it is not being tracked '
                    'by the version logger.')
                version_logger.write_default(key)
            print ('Reusing a cached file.')
            version_logger.update_cached(key)
        # Retrieve a file
        else:
            print ('Retrieving the file...')
            version = self.retrieve_file(temp_filename)
            version_logger.write_retrieved(key, version)


    def retrieve_file(
        self,
        temp_filename: Path,
        *args,
        **kwargs,
    ) -> str:

        # Intended to be overloaded in derived classes
        return 'unknown'