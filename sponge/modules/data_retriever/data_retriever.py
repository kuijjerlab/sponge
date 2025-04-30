### Imports ###
from pathlib import Path
from typing import List

from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .region_retriever import RegionRetriever
from .tfbs_retriever import TFBSRetriever

### Class definition ###
class DataRetriever:
    # Methods
    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):
        """
        Class that retrieves the data required for the running of
        SPONGE, notably the region and TFBS files.

        Parameters
        ----------
        temp_folder : Path
            Folder to save the retrieved file to
        core_config : ConfigManager
            Core configuration of SPONGE
        user_config : ConfigManager
            User-provided configuration of SPONGE
        version_logger : VersionLogger
            Version logger to keep track of the retrieved files
        """

        # JASPAR bigbed file (if appropriate)
        self.tfbs = TFBSRetriever(
            temp_folder,
            core_config['url']['motif']['full'],
            user_config['motif'],
            user_config['genome_assembly'],
            user_config['on_the_fly_processing'],
        )
        # Regions of interest (promoters by default)
        # along with their mapping to genes
        self.regions = RegionRetriever(
            temp_folder,
            core_config['url']['region']['xml'],
            core_config['url']['region']['rest'],
            core_config['url']['chrom_mapping'],
            user_config['region'],
            user_config['genome_assembly'],
            core_config['default_mapping'],
            core_config['default_chromosomes'],
        )
        # Register the classes for version logging
        version_logger.register_class(self.tfbs)
        version_logger.register_class(self.regions)


    def retrieve_data(
        self,
    ) -> None:
        """
        Retrieve the TFBS and region files based on the information
        provided during initialisation.
        """

        self.tfbs.retrieve_file()
        self.regions.retrieve_file()


    def get_tfbs_path(
        self,
    ) -> Path:
        """
        Get the actual path to the TFBS file, can be None if it wasn't
        retrieved.

        Returns
        -------
        Path
            Path to the TFBS file
        """

        return self.tfbs.actual_path


    def get_regions_path(
        self,
    ) -> Path:
        """
        Get the actual path to the region file, can be None if it wasn't
        retrieved.

        Returns
        -------
        Path
            Path to the region file
        """

        return self.regions.actual_path