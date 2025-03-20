### Imports ###
import bioframe
import time

import pandas as pd

from pathlib import Path
from typing import Iterable, Optional

from sponge.config_manager import ConfigManager
from sponge.modules.utils import iterate_chromosomes, iterate_motifs
from sponge.modules.version_logger import VersionLogger

### Class definition ###
class MatchFilter:
    # Functions
    def __init__(
        self,
        core_config: ConfigManager,
        user_config: ConfigManager,
        version_logger: VersionLogger,
        tfbs_path: Optional[Path],
        regions_path: Path,
    ):

        self.url = core_config['url']['motif']['by_tf']
        self.default_chromosomes = core_config['default_chromosomes']
        self.settings = user_config['filter']
        self.assembly = user_config['genome_assembly']
        self.jaspar_release = user_config['motif']['jaspar_release']
        self.version_logger = version_logger
        self.tfbs_path = tfbs_path
        self.regions_path = regions_path


    def filter_matches(
        self,
        tf_names: Iterable[str],
        matrix_ids: Iterable[str],
        score_threshold: Optional[float] = None,
        chromosomes: Optional[Iterable[str]] = None,
        n_processes: Optional[int] = None,
    ) -> None:
        """
        Filters all the binding sites in the JASPAR bigbed file to
        select only the ones in the promoter regions of genes on given
        chromosomes, subject to a score threshold. Stores the result
        internally.

        Parameters
        ----------
        score_threshold : Optional[float], optional
            Minimal score of a match for it to be included in the
            prior or None to follow the option from the initialisation,
            by default None
        chromosomes : Optional[Iterable[str]], optional
            Which chromosomes to get the promoters from or None to
            follow the option from the initialisation, by default None
        n_processes : Optional[int], optional
            Number of processes to run in parallel or None to
            follow the option from the initialisation, by default None
        """

        if n_processes is None:
            n_processes = self.settings['n_processes']
        if chromosomes is None:
            chromosomes = self.default_chromosomes
        if score_threshold is None:
            score_threshold = self.settings['score_threshold']

        print ()
        print ('Loading the regions file...')
        df_full = bioframe.read_table(self.regions_path, header=0)
        df_full.set_index('Transcript stable ID', inplace=True)

        start_time = time.time()
        if self.tfbs_path is None:
            results_list = iterate_motifs(df_full, tf_names,
                matrix_ids, chromosomes, self.jaspar_release,
                self.assembly, n_processes, score_threshold)
            self.version_logger.write_retrieved('tfbs_file',
                self.jaspar_release)
        else:
            results_list = iterate_chromosomes(self.tfbs_path, df_full,
                matrix_ids, chromosomes, n_processes, score_threshold)

        elapsed = time.time() - start_time

        print ()
        print (f'Total time: {elapsed // 60:n} m {elapsed % 60:.2f} s')

        # Save the final results, ignoring the index makes this fast
        # The index is irrelevant
        self.all_edges = pd.concat(results_list, ignore_index=True)