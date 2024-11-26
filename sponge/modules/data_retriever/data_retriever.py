### Imports ###
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.version_logger import VersionLogger

from .motif_retriever import MotifRetriever

### Class definition ###
class DataRetriever:
    
    def __init__(
        self,
        temp_folder: Path,
        config: ConfigReader,
        version_logger: VersionLogger,
    ):
        
        pass