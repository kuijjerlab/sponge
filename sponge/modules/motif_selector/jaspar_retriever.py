### Imports ###
from Bio.motifs.jaspar import Motif
from collections import defaultdict
from pyjaspar import jaspardb
from typing import List, Mapping

from sponge.modules.utils import calculate_ic

### Class definition ###
class JasparRetriever:
    # Methods
    def __init__(
        self,
        motif_settings: dict,
    ):
        """
        Class that retrieves TF motifs from JASPAR.

        Parameters
        ----------
        motif_settings : dict
            Motif relevant settings from the user configuration,
            specifically containing jaspar_release, drop_heterodimers,
            unique_motifs, tf_names and matrix_ids
        """

        self.motif_settings = motif_settings
        # To be overwritten once retrieved
        self.motifs = None
        self.tf_to_motif = None
        # Overwritten by registering with a VersionLogger instance
        self.version_logger = None


    def retrieve_tfs(
        self,
    ) -> None:
        """
        Retrieve the TF motifs of interest based on the settings
        provided during initialisation.
        """

        print ('\n--- Retrieving transcription factor motifs ---')

        jdb_obj = jaspardb(self.motif_settings['jaspar_release'])
        self.write_retrieved('jaspar_motifs',
            self.motif_settings['jaspar_release'])

        drop_heterodimers = self.motif_settings['drop_heterodimers']
        unique_motifs = self.motif_settings['unique_motifs']
        tf_names = self.motif_settings['tf_names']
        matrix_ids = self.motif_settings['matrix_ids']

        if (matrix_ids is not None and len(matrix_ids) > 0 and
            tf_names is not None and len(tf_names) > 0):
            print ('Both motif IDs and TF names have been specified, will '
                'filter on both (intersection)')

        # Latest vertebrate motifs, filter by matrix IDs if any
        motifs = jdb_obj.fetch_motifs(collection='CORE',
            tax_group='vertebrates', matrix_id=matrix_ids)
        # Filter also by TF names if any
        if tf_names is not None and len(tf_names) > 0:
            tf_name_set = set(tf_names)
            motifs_filt = [motif for motif in motifs
                if motif.name in tf_name_set]
        else:
            motifs_filt = motifs
        print ('Retrieved motifs:', len(motifs_filt))

        tf_to_motif = defaultdict(dict)
        for motif in motifs_filt:
            tf_to_motif[motif.name][motif.matrix_id] = calculate_ic(motif)
        self.tf_to_motif = tf_to_motif
        if unique_motifs:
            # Keep only one motif per TF
            motifs_unique = [motif for motif in motifs_filt if
                (tf_to_motif[motif.name][motif.matrix_id] ==
                max(tf_to_motif[motif.name].values()))]
            print ('Unique motifs:', len(motifs_unique))
            motifs_filt = motifs_unique

        if drop_heterodimers:
            motifs_nohd = [motif for motif in motifs_filt
                if '::' not in motif.name]
            print ('Motifs without heterodimers:', len(motifs_nohd))
            self.motifs = motifs_nohd
        else:
            self.motifs = motifs_filt


    def get_motifs(
        self,
    ) -> List[Motif]:
        """
        Get the list of selected TF motifs, can be None if it wasn't
        retrieved.

        Returns
        -------
        List[Motif]
            List of selected TF motifs
        """

        return self.motifs


    def get_tf_to_motif(
        self,
    ) -> Mapping[str, dict]:
        """
        Get the mapping of TF names to their motif IDs and information
        contents, can be None if it wasn't retrieved.

        Returns
        -------
        Mapping[str, dict]
            Mapping of TF names
        """

        return self.tf_to_motif

    # Placeholder functions to be replaced with VersionLogger if any
    def write_retrieved(
        self,
        *args,
        **kwargs,
    ) -> None:

        pass