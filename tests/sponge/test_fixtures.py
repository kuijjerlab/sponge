### Imports ###
import bioframe
import os
import pytest

import pandas as pd

from Bio.motifs.jaspar import Motif
from pyjaspar import jaspardb

from sponge.config_manager import ConfigManager

### Fixtures ###
# Core config fixture
@pytest.fixture
def core_config():
    return ConfigManager()

# A motif without any information
@pytest.fixture
def no_info_motif():
    no_info_row = [0.25] * 4
    no_info_counts = [no_info_row] * 6
    no_info_pwm = pd.DataFrame(no_info_counts, columns=['A', 'C', 'G', 'T'])
    no_info_motif = Motif(matrix_id='XXX', name='XXX', counts=no_info_pwm)

    yield no_info_motif

# A motif with perfect information
@pytest.fixture
def all_A_motif():
    all_A_row = [1] + [0] * 3
    all_A_counts = [all_A_row] * 6
    all_A_pwm = pd.DataFrame(all_A_counts, columns=['A', 'C', 'G', 'T'])
    all_A_motif = Motif(matrix_id='XXX', name='XXX', counts=all_A_pwm)

    yield all_A_motif

# A real motif for SOX2
@pytest.fixture
def SOX2_motif():
    jdb_obj = jaspardb(release='JASPAR2024')
    SOX2_motif = jdb_obj.fetch_motif_by_id('MA0143.1')

    yield SOX2_motif

# A subset of promoters on chromosome 19
@pytest.fixture
def chr19_promoters():
    path_to_file = os.path.join('tests', 'sponge', 'chr19_subset.tsv')
    df = bioframe.read_table(path_to_file, header=0)
    df.set_index('Transcript stable ID', inplace=True)

    yield df

# Part of the FOXF2 track for chromosome 19
@pytest.fixture
def foxf2_chr19():
    path_to_file = os.path.join('tests', 'sponge', 'foxf2_chr19_subset.tsv')
    df = pd.read_csv(path_to_file, sep='\t')

    yield df