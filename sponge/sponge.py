### Imports ###
from pathlib import Path

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class Sponge:

    def __init__(
        self,
        temp_folder: Path = '.sponge_temp/',
        config_file: Path = 'user_config.yaml',
    ):

        self.temp_folder = temp_folder
        # Load the file with internal module inputs
        self.core_config = ConfigManager()
        # Load the user-provided config file (or use defaults if it doesn't
        # exist)
        self.user_config = ConfigManager(config_file, temp_folder)
        self.version_logger = VersionLogger(temp_folder)

        self.input_data = DataRetriever(temp_folder, self.core_config,
            self.user_config, self.version_logger)
        
        # TODO: Use the retrieved data
        # Motif selection

        # Match filtering

        # PPI retrieval

        # File output


    def retrieve_data(
        self,
        **kwargs,
    ):
        
        self.input_data = DataRetriever(self.temp_folder, self.core_config,
            self.user_config, self.version_logger)