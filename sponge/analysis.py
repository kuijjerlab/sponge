### Imports ###
import os

import pandas as pd

from typing import Union

from sklearn.metrics import classification_report, confusion_matrix

FILE_LIKE = Union[str, bytes, os.PathLike]
MOTIF_COLS = ['tf', 'gene', 'edge']

### Functions ###
def load_prior(
    path: FILE_LIKE,
) -> pd.DataFrame:
    
    return pd.read_csv(path, sep='\t', header=None, names=MOTIF_COLS)


def describe_prior(
    prior: pd.DataFrame,
) -> None:
    
    n_tfs = prior['tf'].nunique()
    n_genes = prior['gene'].nunique()

    print ('Number of unique TFs:', n_tfs)
    print ('Number of unique genes:', n_genes)
    print ('Number of edges:', len(prior))
    print (f'Network density: {100 * len(prior) / (n_tfs * n_genes):.2f} %')


def compare_priors(
    prior_1: pd.DataFrame,
    prior_2: pd.DataFrame,
) -> None:

    print ('Statistics for the first prior:')
    describe_prior(prior_1)
    print ()
    print ('Statistics for the second prior:')
    describe_prior(prior_2)

    common_tfs = set(prior_1['tf'].unique()).intersection(
        prior_2['tf'].unique())
    common_genes = set(prior_1['gene'].unique()).intersection(
        prior_2['gene'].unique())
    print ()
    print ('Number of common TFs:', len(common_tfs))
    print ('Number of common genes:', len(common_genes))

    common_index = pd.MultiIndex.from_product([sorted(common_tfs), 
        sorted(common_genes)])
    prior_1_mod = prior_1.set_index(['tf', 'gene']).reindex(
        common_index, fill_value=0)
    prior_2_mod = prior_2.set_index(['tf', 'gene']).reindex(
        common_index, fill_value=0)
    comp_df = prior_1_mod.join(prior_2_mod, lsuffix='_1', rsuffix='_2')

    print ()
    print ('Network density in common TF/genes for the first prior:',
        f'{100 * comp_df["edge_1"].mean():.2f} %')
    print ('Network density in common TF/genes for the second prior:',
        f'{100 * comp_df["edge_2"].mean():.2f} %')
    print ()
    print (classification_report(comp_df['edge_1'], comp_df['edge_2']))
    print (confusion_matrix(comp_df['edge_1'], comp_df['edge_2']))