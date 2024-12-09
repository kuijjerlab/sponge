### Imports ###
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class Sponge:

    def __init__(
        self,
        temp_folder: Path = '.sponge_temp/',
        config_file: Path = 'user_config.yaml',
    ):

        # Load the file with internal module inputs
        self.core_config = ConfigReader()
        # Load the user-provided config file (or use defaults if it doesn't
        # exist)
        self.user_config = ConfigReader(config_file, temp_folder)
        self.version_logger = VersionLogger(temp_folder)

        self.input_data = DataRetriever(temp_folder, self.core_config,
            self.user_config, self.version_logger)