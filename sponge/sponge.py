### Imports ###
from pathlib import Path

from sponge.config_reader import ConfigReader
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class Sponge:

    def __init__(
        self,
        temp_folder: Path = '.sponge_temp/',
    ):

        self.config = ConfigReader().config
        self.version_logger = VersionLogger()

        self.input_data = DataRetriever(self.logger)