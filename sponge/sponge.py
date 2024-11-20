### Imports ###
from sponge.config_reader import ConfigReader
from sponge.modules.data_retriever import DataRetriever
from sponge.modules.logger import Logger

### Class definition ###
class Sponge:

    def __init__(
        self,
    ):

        self.config = ConfigReader().config
        self.logger = Logger()

        self.input_data = DataRetriever(self.logger)