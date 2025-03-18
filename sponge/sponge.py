### Imports ###
from pathlib import Path

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.match_filter import MatchFilter
from sponge.modules.motif_selector import MotifSelector
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class Sponge:
    # Functions
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

        data = DataRetriever(self.temp_folder, self.core_config,
            self.user_config, self.version_logger)

        self.tfbs_path = data.tfbs.actual_path
        self.regions_path = data.regions.actual_path
        self.chromosomes = data.regions.settings['chromosomes']


    def select_motifs(
        self,
        **kwargs,
    ):

        motifs = MotifSelector(self.core_config, self.user_config,
            self.version_logger)
        motifs.select_tfs()
        motifs.find_homologs()

        self.motifs = motifs.motifs
        self.animal_to_human = motifs.animal_to_human
        self.matrix_ids = motifs.matrix_ids
        self.tf_names = motifs.tf_names


    def filter_tfbs(
        self,
        **kwargs,
    ):

        match_filter = MatchFilter(self.core_config, self.user_config,
            self.version_logger, self.tfbs_path, self.regions_path)
        match_filter.filter_matches(self.tf_names, self.matrix_ids,
            chromosomes=self.chromosomes)

        self.all_edges = match_filter.all_edges


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