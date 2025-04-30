### Imports ###
from bioframe import assembly_info
from functools import wraps
from pathlib import Path
from typing import Union

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.file_writer import FileWriter
from sponge.modules.match_aggregator import MatchAggregator
from sponge.modules.match_filter import MatchFilter
from sponge.modules.motif_selector import MotifSelector
from sponge.modules.ppi_retriever import PPIRetriever
from sponge.modules.utils import process_jaspar_version
from sponge.modules.version_logger import VersionLogger

# General TODOs:
# TODO: For most class definitions, refactor out passing the entire config
# (just the relevant parts)
# TODO: Potentially refactor out VersionLogger too
# Ideally a wrapper would be implemented that would log the version
# Requires the classes to store the version information after retrieval
# Would this mess up the accuracy of the timestamp?

### Decorators ###
def allow_config_update(
    function,
) -> None:

    # sig = utilipy.utils.inspect.fuller_signature(function)
    # param = inspect.Parameter('kw1', inspect._KEYWORD_ONLY, default='has a value', annotation='added kwarg')
    # sig = sig.insert_parameter(sig.index_var_keyword, param)

    # @utilipy.utils.functools.wraps(function, signature=sig)

    @wraps(function)
    def wrapper(
        self,
        user_config_update: dict = {}
    ) -> None:

        self.user_config.deep_update(user_config_update)
        self.fill_default_values()
        function(self)

    # TODO: Find a way to update signature too
    if wrapper.__doc__ is None:
        wrapper.__doc__ = ''
    wrapper.__doc__ += ("\nbbb")

    return wrapper

### Class definition ###
class Sponge:
    # Methods
    def __init__(
        self,
        temp_folder: Path = '.sponge_temp/',
        config: Union[Path, dict] = 'user_config.yaml',
        config_update: dict = {},
    ):

        self.temp_folder = temp_folder
        # Load the file with internal module inputs
        self.core_config = ConfigManager()
        # Load the user-provided config file
        # (or use defaults if it doesn't exist)
        self.user_config = ConfigManager(config, temp_folder)
        self.user_config.deep_update(config_update)
        # Fill in the default values in the user config
        self.fill_default_values()
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


    def fill_default_values(
        self,
    ) -> None:

        # Fill in chromosomes
        if self.user_config['region']['use_all_chromosomes']:
            all_names = list(assembly_info(self.user_config['genome_assembly'],
                roles='all').seqinfo['name'])
            self.user_config['region']['chromosomes'] = all_names
        elif self.user_config['region']['chromosomes'] is None:
            d_chromosomes = self.core_config['default_chromosomes']
            self.user_config['region']['chromosomes'] = d_chromosomes
        # Fill in JASPAR release
        if self.user_config['motif']['jaspar_release'] is None:
            newest_release = process_jaspar_version(None)
            self.user_config['motif']['jaspar_release'] = newest_release


    def retrieve_data(
        self,
    ) -> None:

        data = DataRetriever(
            self.temp_folder,
            self.core_config,
            self.user_config,
            self.version_logger,
        )
        data.retrieve_data()

        self.tfbs_path = data.get_tfbs_path()
        self.regions_path = data.get_regions_path()


    @allow_config_update
    def select_motifs(
        self,
    ) -> None:
        """
        aaa
        """

        motifs = MotifSelector(
            self.core_config,
            self.user_config,
            self.version_logger,
        )
        motifs.select_tfs()
        motifs.find_homologs()

        self.homolog_mapping = motifs.get_homolog_mapping()
        self.matrix_ids = motifs.get_matrix_ids()
        self.tf_names = motifs.get_tf_names()


    @allow_config_update
    def filter_tfbs(
        self,
    ) -> None:

        match_filter = MatchFilter(self.core_config, self.user_config,
            self.version_logger, self.tfbs_path, self.regions_path)
        match_filter.filter_matches(
            self.tf_names,
            self.matrix_ids,
            self.user_config['filter']['score_threshold'],
            self.user_config['region']['chromosomes'],
            self.user_config['filter']['n_processes'],
        )

        self.all_edges = match_filter.all_edges
        self.all_edges['weight'] = self.all_edges['score'] / 100
        self.all_edges['edge'] = 1


    @allow_config_update
    def retrieve_ppi(
        self,
    ) -> None:

        ppi = PPIRetriever(self.core_config, self.user_config,
            self.version_logger)
        filtered_tfs = self.all_edges['TFName'].unique()
        mapped_tfs = [self.homolog_mapping[x] if x in self.homolog_mapping
            else x for x in filtered_tfs]
        ppi.retrieve_ppi(mapped_tfs)

        self.ppi_frame = ppi.ppi_frame
        self.ppi_frame['edge'] = 1


    @allow_config_update
    def write_output_files(
        self,
    ) -> None:

        aggregator = MatchAggregator(self.all_edges, self.regions_path,
            self.homolog_mapping)
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