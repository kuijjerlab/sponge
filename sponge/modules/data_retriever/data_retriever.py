### Imports ###
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.version_logger import VersionLogger

from .motif_retriever import MotifRetriever
from .region_retriever import RegionRetriever

### Class definition ###
class DataRetriever:

    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigReader,
        user_config: ConfigReader,
        version_logger: VersionLogger,
    ):

        # Retrieve the following:
        # JASPAR bigbed file if appropriate
        # Regions of interest (promoters) along with their mapping to genes
        self.tfbs = MotifRetriever(temp_folder, core_config, user_config,
            version_logger)
        self.regions = RegionRetriever(temp_folder, core_config, user_config,
            version_logger)