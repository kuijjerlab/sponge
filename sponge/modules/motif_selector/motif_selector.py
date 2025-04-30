### Imports ###
from typing import List, Mapping

from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger

from .homology_retriever import HomologyRetriever
from .jaspar_retriever import JasparRetriever

### Class definition ###
class MotifSelector:
    # Methods
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):
        """
        Class that retrieves the TF motifs from JASPAR and homologs
        from NCBI.

        Parameters
        ----------
        core_config : ConfigManager
            Core configuration of SPONGE
        user_config : ConfigManager
            User-provided configuration of SPONGE
        version_logger : VersionLogger
            Version logger to keep track of the retrieved files
        """

        self.core_config = core_config
        self.user_config = user_config
        self.version_logger = version_logger
        # To be overwritten once retrieved
        self.homolog_mapping = None
        self.matrix_ids = None
        self.tf_names = None


    def select_tfs(
        self,
    ) -> None:
        """
        Selects an initial selection of TFs based on the settings
        provided during initialisation.
        """

        jaspar = JasparRetriever(self.user_config['motif'])
        self.version_logger.register_class(jaspar)
        jaspar.retrieve_tfs()

        self.motifs = jaspar.get_motifs()
        self.tf_to_motif = jaspar.get_tf_to_motif()


    def find_homologs(
        self,
    ) -> None:
        """
        Finds the homologs of non species matching TFs in the species
        of interest using the settings provided during initialisation.
        """

        homologs = HomologyRetriever(
            self.user_config['motif']['unique_motifs'],
            self.core_config['url']['protein'],
            self.core_config['url']['homology'],
        )
        self.version_logger.register_class(homologs)
        homologs.find_homologs(self.motifs, self.tf_to_motif)

        self.homolog_mapping = homologs.get_homolog_mapping()
        self.matrix_ids = homologs.get_matrix_ids()
        self.tf_names = homologs.get_tf_names()


    def get_homolog_mapping(
        self,
    ) -> Mapping[str, str]:
        """
        Get the mapping of TF homolog names to the names in the species
        of interest, can be None if it wasn't retrieved.

        Returns
        -------
        Mapping[str, str]
            Mapping of TF homolog names
        """

        return self.homolog_mapping


    def get_matrix_ids(
        self,
    ) -> List[str]:
        """
        Get the list of TF matrix IDs, can be None if it wasn't
        retrieved.

        Returns
        -------
        List[str]
            List of TF matrix IDs
        """

        return self.matrix_ids


    def get_tf_names(
        self,
    ) -> List[str]:
        """
        Get the list of TF names, can be None if it wasn't retrieved.

        Returns
        -------
        List[str]
            List of TF names
        """

        return self.tf_names