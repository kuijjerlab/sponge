### Imports ###
from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .homology_retriever import HomologyRetriever
from .jaspar_retriever import JasparRetriever
from .protein_id_mapper import ProteinIDMapper

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

        mapper = ProteinIDMapper(self.core_config)
        # Retrieve mapping of Uniprot to GeneID
        all_ids = set()
        for motif in self.motifs:
            for id in motif.acc:
                all_ids.add(id)
        mapping = mapper.get_uniprot_mapping('UniProtKB_AC-ID', 'GeneID',
            list(all_ids))
        homologs = HomologyRetriever(self.core_config, self.user_config,
            self.version_logger)
        homologs.find_homologs(self.motifs, self.tf_to_motif, mapping)

        self.animal_to_human = homologs.animal_to_human
        self.matrix_ids = homologs.matrix_ids
        self.tf_names = homologs.tf_names
        