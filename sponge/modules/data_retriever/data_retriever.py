### Imports ###
from pathlib import Path

from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .region_retriever import RegionRetriever
from .tfbs_retriever import TFBSRetriever

### Class definition ###
class DataRetriever:
    # Functions
    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):

        # Retrieve the following:
        # JASPAR bigbed file if appropriate
        self.tfbs = TFBSRetriever(temp_folder, core_config, user_config,
            version_logger)
        # Regions of interest (promoters) along with their mapping to genes
        self.regions = RegionRetriever(temp_folder, core_config, user_config,
            version_logger)