import pandas as pd

from Bio.motifs.jaspar import Motif

from pyjaspar import jaspardb

# A motif without any information
no_info_row = [0.25] * 4
no_info_counts = [no_info_row] * 6
no_info_pwm = pd.DataFrame(no_info_counts, columns=['A', 'C', 'G', 'T'])
no_info_motif = Motif(matrix_id='XXX', name='XXX', counts=no_info_pwm)

# A motif with perfect information
all_A_row = [1] + [0] * 3
all_A_counts = [all_A_row] * 6
all_A_pwm = pd.DataFrame(all_A_counts, columns=['A', 'C', 'G', 'T'])
all_A_motif = Motif(matrix_id='XXX', name='XXX', counts=all_A_pwm)

# A real motif for SOX2
jdb_obj = jaspardb(release='JASPAR2024')
SOX2_motif = jdb_obj.fetch_motif_by_id('MA0143.1')