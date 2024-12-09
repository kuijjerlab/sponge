### Imports ###
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class RegionRetriever:
    
    def __init__(
        self,
        temp_folder: Path,
        core_config: ConfigReader,
        user_config: ConfigReader,
        version_logger: VersionLogger,
    ):

        pass