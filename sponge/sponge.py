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
        # Retrieve necessary files if required
        self.retrieve_data()
        # Run the default workflow if selected
        if not self.user_config.is_false('default_workflow'):
            self.select_motifs()
            self.filter_tfbs()
            self.retrieve_ppi()
            self.write_output_files()
        # Otherwise, let the user call the functions individually


    def retrieve_data(
        self,
        **kwargs,
    ):
        
        self.input_data = DataRetriever(self.temp_folder, self.core_config,
            self.user_config, self.version_logger)
        
    
    def select_motifs(
        self,
        **kwargs,
    ):
        
        pass


    def filter_tfbs(
        self,
        **kwargs,
    ):
        
        pass


    def retrieve_ppi(
        self,
        **kwargs,
    ):
        
        pass


    def write_output_files(
        self,
        **kwargs,
    ):
        
        pass