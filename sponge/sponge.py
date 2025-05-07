### Imports ###
import inspect

from bioframe import assembly_info
from functools import wraps
from pathlib import Path
from typing import Callable, Union

from sponge.config_manager import ConfigManager
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.file_writer import FileWriter
from sponge.modules.match_aggregator import MatchAggregator
from sponge.modules.match_filter import MatchFilter
from sponge.modules.motif_selector import MotifSelector
from sponge.modules.ppi_retriever import PPIRetriever
from sponge.modules.utils import process_jaspar_version
from sponge.modules.version_logger import VersionLogger

### Decorators ###
def allow_config_update(
    function: Callable,
) -> Callable:
    """
    Adds the ability for a method to update the user configuration
    through an optional argument user_config_update or directly through
    keywords.

    Parameters
    ----------
    function : Callable
        Method to be modified

    Returns
    -------
    Callable
        Modified version of the method
    """

    @wraps(function)
    def wrapper(
        self,
        user_config_update: dict = {},
        **kwargs,
    ) -> None:

        self.user_config.deep_update(user_config_update)
        self.user_config.deep_update(kwargs)
        self.fill_default_values()
        function(self)

    # Update function signature
    orig_signature = inspect.signature(function)
    new_params = orig_signature.parameters.copy()
    new_params['user_config_update'] = inspect.Parameter('user_config_update',
        inspect.Parameter.POSITIONAL_OR_KEYWORD, default={})
    new_params['**kwargs'] = inspect.Parameter('kwargs',
        inspect.Parameter.VAR_KEYWORD)
    new_signature = inspect.Signature(parameters=new_params.values())
    wrapper.__signature__ = new_signature
    # Update function docstring
    if wrapper.__doc__ is None:
        wrapper.__doc__ = ''
    wrapper.__doc__ += (
        """
        Parameters
        ----------
        user_config_update : dict, optional
            Dictionary to be used for updating the user configuration
            prior to executing the method, by default {}
        **kwargs
            Keyword arguments will be used to update the user
            configuration in the class instance prior to executing
            the method, they take priority over user_config_update
        """
    )

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
        """
        Class that generates a prior gene regulatory network and a prior
        protein-protein interaction network based on the provided
        settings.

        Parameters
        ----------
        temp_folder : Path, optional
            Folder to save the temporary files to,
            by default '.sponge_temp/'
        config : Union[Path, dict], optional
            Path to a configuration file in .yaml format or a dictionary
            containing the settings, by default 'user_config.yaml'
        config_update : dict, optional
            Dictionary that updates the provided settings, useful for
            changing a small subset of the settings from a config file,
            by default {}
        """

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
        """
        Fills in the values for chromosomes and jaspar_release in the
        stored user configuration if needed.
        """

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
        """
        Retrieves the initial files required for running the SPONGE
        pipeline. Some online services are still accessed later
        in the process.
        """

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
        Selects the TF motifs to be used in the SPONGE pipeline based
        on the settings from the configuration.
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
        """
        Filter the TF binding sites to match the regions of interest
        and score threshold from the configuration.
        """

        match_filter = MatchFilter(
            self.tfbs_path,
            self.regions_path,
            self.user_config['genome_assembly'],
            self.user_config['motif']['jaspar_release'],
            self.core_config['url']['motif']['by_tf'],
        )
        self.version_logger.register_class(match_filter)
        match_filter.filter_matches(
            self.user_config['region']['chromosomes'],
            self.matrix_ids,
            self.tf_names,
            self.user_config['filter']['score_threshold'],
            self.user_config['filter']['n_processes'],
        )

        self.all_edges = match_filter.get_all_edges()
        self.all_edges['weight'] = self.all_edges['score'] / 100
        self.all_edges['edge'] = 1


    @allow_config_update
    def retrieve_ppi(
        self,
    ) -> None:
        """
        Retrieves PPI data for selected TFs from the STRING database.
        """

        ppi = PPIRetriever(
            self.core_config['url']['ppi'],
            self.core_config['url']['protein'],
        )
        self.version_logger.register_class(ppi)
        filtered_tfs = self.all_edges['TFName'].unique()
        mapped_tfs = [self.homolog_mapping[x] if x in self.homolog_mapping
            else x for x in filtered_tfs]
        ppi.retrieve_ppi(
            mapped_tfs,
            self.user_config['ppi']['physical_only'],
        )

        self.ppi_frame = ppi.get_ppi_frame()
        self.ppi_frame['edge'] = 1


    @allow_config_update
    def write_output_files(
        self,
    ) -> None:
        """
        Writes the generated gene regulatory and PPI priors to files.
        """

        aggregator = MatchAggregator(
            self.all_edges,
            self.regions_path,
            self.homolog_mapping,
        )
        aggregator.aggregate_matches(
            self.user_config['motif_output']['use_gene_names'],
            self.user_config['motif_output']['protein_coding_only'],
        )
        edges = aggregator.get_edges()

        writer = FileWriter()

        print ('\n--- Saving the motif prior ---')
        motif_weight = 'edge'
        if self.user_config['motif_output']['weighted']:
            motif_weight = 'weight'
        label = 'Gene stable ID'
        if self.user_config['motif_output']['use_gene_names']:
            label = 'Gene name'
        writer.write_network_file(
            edges,
            ['TFName', label],
            motif_weight,
            self.user_config['motif_output']['file_name'],
        )

        print ('\n--- Saving the PPI prior ---')
        ppi_weight = 'edge'
        if self.user_config['ppi_output']['weighted']:
            ppi_weight = 'score'
        writer.write_network_file(
            self.ppi_frame,
            ['tf1', 'tf2'],
            ppi_weight,
            self.user_config['ppi_output']['file_name'],
        )