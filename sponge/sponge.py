### Imports ###
from pathlib import Path
from typing import Union

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.file_writer import FileWriter
from sponge.modules.match_aggregator import MatchAggregator
from sponge.modules.match_filter import MatchFilter
from sponge.modules.motif_selector import MotifSelector
from sponge.modules.ppi_retriever import PPIRetriever
from sponge.modules.version_logger import VersionLogger

# General TODOs:
# TODO: For most class definitions, refactor out passing the entire config
# (just the relevant parts)
# TODO: Potentially refactor out VersionLogger too
# Ideally a wrapper would be implemented that would log the version
# Requires the classes to store the version information after retrieval
# Would this mess up the accuracy of the timestamp?

### Class definition ###
class Sponge:
    # Methods
    def __init__(
        self,
        temp_folder: Path = '.sponge_temp/',
        config: Union[Path, dict] = 'user_config.yaml',
    ):

        self.temp_folder = temp_folder
        # Load the file with internal module inputs
        self.core_config = ConfigManager()
        # Load the user-provided config file
        # (or use defaults if it doesn't exist)
        self.user_config = ConfigManager(config, temp_folder)
        # TODO: Fill in some default values at this stage?
        self.version_logger = VersionLogger(temp_folder)
        # Retrieve necessary files if required
        self.retrieve_data()
        # Run the default workflow if selected
        if self.user_config['default_workflow']:
            self.select_motifs()
            self.filter_tfbs()
            self.retrieve_ppi()
            self.write_output_files()
        # Otherwise, let the user call the functions individually

    # TODO: Add proper arguments to the functions to potentially overwrite
    # the ones from the config file
    def retrieve_data(
        self,
        **kwargs,
    ):

        data = DataRetriever(self.temp_folder, self.core_config,
            self.user_config, self.version_logger)
        data.retrieve_data()

        self.tfbs_path = data.get_tfbs_path()
        self.regions_path = data.get_regions_path()
        self.chromosomes = data.get_chromosomes()


    def select_motifs(
        self,
        **kwargs,
    ):

        motifs = MotifSelector(self.core_config, self.user_config,
            self.version_logger)
        motifs.select_tfs()
        motifs.find_homologs()

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
        self.all_edges['weight'] = self.all_edges['score'] / 100
        self.all_edges['edge'] = 1


    def retrieve_ppi(
        self,
        **kwargs,
    ):

        ppi = PPIRetriever(self.core_config, self.user_config,
            self.version_logger)
        filtered_tfs = self.all_edges['TFName'].unique()
        humanised_tfs = [self.animal_to_human[x] if x in self.animal_to_human
            else x for x in filtered_tfs]
        ppi.retrieve_ppi(humanised_tfs)

        self.ppi_frame = ppi.ppi_frame
        self.ppi_frame['edge'] = 1


    def write_output_files(
        self,
        **kwargs,
    ):

        aggregator = MatchAggregator(self.all_edges, self.regions_path,
            self.animal_to_human)
        aggregator.aggregate_matches(
            self.user_config['motif_output']['use_gene_names'],
            self.user_config['motif_output']['protein_coding_only'],
        )
        edges = aggregator.edges

        writer = FileWriter()

        print ('\n--- Saving the motif prior ---')
        motif_weight = 'edge'
        if self.user_config['motif_output']['weighted']:
            motif_weight = 'weight'
        label = 'Gene stable ID'
        if self.user_config['motif_output']['use_gene_names']:
            label = 'Gene name'
        writer.write_network_file(edges, ['TFName', label], motif_weight,
            self.user_config['motif_output']['file_name'])

        print ('\n--- Saving the PPI prior ---')
        ppi_weight = 'edge'
        if self.user_config['ppi_output']['weighted']:
            ppi_weight = 'score'
        writer.write_network_file(self.ppi_frame, ['tf1', 'tf2'], ppi_weight,
            self.user_config['ppi_output']['file_name'])