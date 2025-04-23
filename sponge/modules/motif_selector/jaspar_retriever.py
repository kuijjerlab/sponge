### Imports ###
from collections import defaultdict
from pyjaspar import jaspardb
from typing import Optional

from sponge.config_manager import ConfigManager
from sponge.modules.version_logger import VersionLogger
from sponge.modules.utils import calculate_ic

### Class definition ###
class JasparRetriever:
    # Methods
    def __init__(
        self,
        user_config: ConfigManager,
        version_logger: VersionLogger,
    ):

        self._initialise_jaspar(user_config.get_value(
            ['motif', 'jaspar_release']))
        version_logger.write_retrieved('jaspar_motifs', self.jaspar_release)

        self.user_config = user_config
        self.version_logger = version_logger


    def _initialise_jaspar(
        self,
        jaspar_release: Optional[str] = None,
    ) -> None:
        """
        Initialises the JASPAR database object with the best matching
        release.

        Parameters
        ----------
        jaspar_release : Optional[str], optional
            Which JASPAR release to use or None to select the newest,
            by default None
        """

        # Initialise a database object to interact with
        self.jdb_obj = jaspardb()
        if jaspar_release is None:
            # Just keep the current object, log the release version
            from pyjaspar import JASPAR_LATEST_RELEASE
            self.jaspar_release = JASPAR_LATEST_RELEASE
        else:
            jaspar_available = self.jdb_obj.get_releases()
            alt_name = 'JASPAR' + jaspar_release
            if jaspar_release in jaspar_available:
                # Release found as specified
                self.jaspar_release = jaspar_release
                self.jdb_obj = jaspardb(self.jaspar_release)
            elif alt_name in jaspar_available:
                # Try adding JASPAR to the provided release
                # Converts e.g. 2022 to JASPAR2022 (actual release name)
                print (f'Found {alt_name} in available releases, assuming '
                    'this matches your choice')
                self.jaspar_release = alt_name
                self.jdb_obj = jaspardb(self.jaspar_release)
            else:
                error_str = ('The specified JASPAR release '
                    f'({jaspar_release}) is not available.\n'
                    'Available releases: ' +
                    ', '.join(jaspar_available))
                raise ValueError(error_str)


    def retrieve_tfs(
        self,
    ):

        print ('\n--- Retrieving transcription factor motifs ---')

        drop_heterodimers = self.user_config['motif']['drop_heterodimers']
        unique_motifs = self.user_config['motif']['unique_motifs']
        tf_names = None
        if self.user_config.exists(['motif', 'tf_names']):
            tf_names = self.user_config['motif']['tf_names']
        matrix_ids = None
        if self.user_config.exists(['motif', 'matrix_ids']):
            matrix_ids = self.user_config['motif']['matrix_ids']

        if (matrix_ids is not None and len(matrix_ids) > 0 and
            tf_names is not None and len(tf_names) > 0):
            print ('Both motif IDs and TF names have been specified, will '
                'filter on both (intersection)')

        # Latest vertebrate motifs, filter by matrix IDs if any
        motifs = self.jdb_obj.fetch_motifs(collection='CORE',
            tax_group='vertebrates', matrix_id=matrix_ids)
        # Filter also by TF names if any
        if tf_names is not None and len(tf_names) > 0:
            tf_name_set = set(tf_names)
            motifs_filt = [i for i in motifs if i.name in tf_name_set]
        else:
            motifs_filt = motifs
        print ('Retrieved motifs:', len(motifs_filt))


        tf_to_motif = defaultdict(dict)
        for i in motifs_filt:
            tf_to_motif[i.name][i.matrix_id] = calculate_ic(i)
        self.tf_to_motif = tf_to_motif
        # Keep only one motif per TF
        if unique_motifs:
            motifs_unique = [i for i in motifs_filt if
                (tf_to_motif[i.name][i.matrix_id] ==
                max(tf_to_motif[i.name].values()))]
            print ('Unique motifs:', len(motifs_unique))
            motifs_filt = motifs_unique

        # Drop heterodimers
        if drop_heterodimers:
            motifs_nohd = [i for i in motifs_filt if '::' not in i.name]
            print ('Motifs without heterodimers:', len(motifs_nohd))
            self.motifs = motifs_nohd
        else:
            self.motifs = motifs_filt