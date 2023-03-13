### Imports ###
import pandas as pd
import bioframe

from math import log2

### Data ###
# TODO: Add URLs
file_df = pd.DataFrame(
    {'description': ['homologene', 'promoter', 'jaspar_bigbed'],
     'name': ['homologene.tsv', 'promoters.bed', 'JASPAR.bb'],
     'url': ['https://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data',
             '',
             ''],
     'mode': ['w', 'w', 'wb']}
).set_index('description')

### Functions ###
def filter_edges(bb_ref, bed_df, motif_list, chrom, start_ind, final_ind, score_threshold=400):
    df = pd.DataFrame()
    for transcript in bed_df.index[start_ind:final_ind]:
        start,end = bed_df.loc[transcript][['start', 'end']]
        motifs = bioframe.read_bigbed(bb_ref, chrom, start=start, end=end)
        max_scores = motifs[motifs['TFName'].isin(motif_list)].groupby('TFName')[['score', 'name']].max()
        max_scores = max_scores[max_scores['score'] >= score_threshold]
        max_scores.reset_index(inplace=True)
        max_scores['transcript'] = transcript
        df = pd.concat((df, max_scores), ignore_index=True, copy=False)
    return df


def filter_edges_helper(input_tuple):
    return filter_edges(*input_tuple)


def plogp(x):
    """
    Returns x*log2(x) for a number, handles the 0 case properly.
    """
    if x == 0:
        return 0
    else:
        return x*log2(x)


def calculate_ic(motif):
    """
    Calculates the information content for a given motif, assumes equal ACGT distribution.
    """
    df = pd.DataFrame(motif.pwm)
    df['IC'] = df.apply(lambda x: 2 + sum([plogp(x[y]) for y in ['A', 'C', 'G', 'T']]), axis=1)
    
    return df['IC'].sum()


def adjust_name(i):
    return i[:-2] + i[-2:].lower()