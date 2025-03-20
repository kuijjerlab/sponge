### Imports ###
from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .homology_retriever import HomologyRetriever
from .jaspar_retriever import JasparRetriever

### Class definition ###
class MotifSelector:
    # Functions
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):

        self.core_config = core_config
        self.user_config = user_config
        self.version_logger = version_logger


    def select_tfs(
        self,
    ):

        jaspar = JasparRetriever(self.user_config, self.version_logger)
        jaspar.retrieve_tfs()

        self.motifs = jaspar.motifs
        self.tf_to_motif = jaspar.tf_to_motif


    def find_homologs(
        self,
    ):

        homologs = HomologyRetriever(self.core_config, self.user_config,
            self.version_logger)
        homologs.find_homologs(self.motifs, self.tf_to_motif)

        self.animal_to_human = homologs.animal_to_human
        self.matrix_ids = homologs.matrix_ids
        self.tf_names = homologs.tf_names
        