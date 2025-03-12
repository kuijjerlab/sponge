### Imports ###
from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .homology_retriever import HomologyRetriever
from .jaspar_retriever import JasparRetriever
from .protein_retriever import ProteinRetriever

### Class definition ###
class MotifSelector:
    # Functions
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):

        pass


    def select_tfs(
        self,
    ):

        pass


    def find_homologs(
        self,
    ):

        pass