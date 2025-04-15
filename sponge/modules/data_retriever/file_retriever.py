### Imports ###
import os

from pathlib import Path
from typing import Callable, Optional

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

        # print ()
        # print (f'--- Attempting to locate {key} ---')

        # # Check existence of a user-provided file
        # if path_to_file is not None:
        #     print (f'Using a user-provided file: {path_to_file}')
        #     if not os.path.exists(path_to_file):
        #         raise FileNotFoundError('Could not locate file: '
        #             f'{path_to_file}')
        #     version_logger.write_provided(key)
        #     self.actual_path = path_to_file
        # # Check for a cached file
        # elif os.path.exists(temp_filename):
        #     if key not in version_logger:
        #         print ('A cached file is present but it is not being tracked '
        #             'by the version logger.')
        #         version_logger.write_default(key)
        #     print ('Reusing a cached file.')
        #     version_logger.update_cached(key)
        #     self.actual_path = temp_filename
        # # Retrieve a file
        # else:
        #     print ('Retrieving the file...')
        #     version = self.retrieve_file(temp_filename)
        #     version_logger.write_retrieved(key, version)
        #     self.actual_path = temp_filename

        self.key = key
        self.temp_filename = temp_filename
        self.version_logger = version_logger
        self.path_to_file = path_to_file


    def retrieve_file(
        self,
        retrieve_function: Callable,
    ) -> str:

        print ()
        print (f'--- Attempting to locate {self.key} ---')

        # Check existence of a user-provided file
        if self.path_to_file is not None:
            print (f'Using a user-provided file: {self.path_to_file}')
            if not os.path.exists(self.path_to_file):
                raise FileNotFoundError('Could not locate file: '
                    f'{self.path_to_file}')
            self.version_logger.write_provided(self.key)
            self.actual_path = self.path_to_file
        # Check for a cached file
        elif os.path.exists(self.temp_filename):
            if self.key not in self.version_logger:
                print ('A cached file is present but it is not being tracked '
                    'by the version logger.')
                self.version_logger.write_default(self.key)
            print ('Reusing a cached file.')
            self.version_logger.update_cached(self.key)
            self.actual_path = self.temp_filename
        # Retrieve a file
        else:
            print ('Retrieving the file...')
            version = retrieve_function()
            self.version_logger.write_retrieved(self.key, version)
            self.actual_path = self.temp_filename