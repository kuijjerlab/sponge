### Imports ###
from collections import defaultdict
from pyjaspar import jaspardb

from sponge.modules.utils import calculate_ic

### Class definition ###
class JasparRetriever:
    # Functions
    def __init__(
        self,
    ):

        pass


    def retrieve_tfs(
        self,
    ):

        if drop_heterodimers is None:
            drop_heterodimers = self.drop_heterodimers
        if tf_names is None:
            tf_names = self.prov_tf_names
        if matrix_ids is None:
            matrix_ids = self.prov_matrix_ids

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

        # Keep only one motif per TF
        # Consider dropping this requirement maybe
        tf_to_motif = defaultdict(dict)
        for i in motifs_filt:
            tf_to_motif[i.name][i.matrix_id] = calculate_ic(i)
        self.tf_to_motif = tf_to_motif
        motifs_unique = [i for i in motifs_filt if
            (tf_to_motif[i.name][i.matrix_id] ==
            max(tf_to_motif[i.name].values()))]
        print ('Unique motifs:', len(motifs_unique))

        # Drop heterodimers
        if drop_heterodimers:
            motifs_nohd = [i for i in motifs_unique if '::' not in i.name]
            print ('Motifs without heterodimers:', len(motifs_nohd))
            self.motifs = motifs_nohd
        else:
            self.motifs = motifs_unique