### Imports ###
import bioframe
import os

import pandas as pd

from typing import Union, Iterable, Tuple

FILE_LIKE = Union[str, bytes, os.PathLike]
FILTER_INPUT = Tuple[str, pd.DataFrame, Iterable[str], str, int, int, float]

### Functions ###
def filter_edges(
    bb_ref: FILE_LIKE, 
    bed_df: pd.DataFrame, 
    motif_list: Iterable[str], 
    chrom: str, 
    start_ind: int, 
    final_ind: int, 
    score_threshold: float = 400
) -> pd.DataFrame:
    """
    Filters possible binding site matches for the provided transcription
    factors from a bigbed file. This is done for a number of continuous
    regions of a single chromosome subject to a score threshold.

    Parameters
    ----------
    bb_ref : FILE_LIKE
        The path to a bigbed file that stores all possible matches
    bed_df : pd.DataFrame
        A pandas DataFrame containing the regions of interest in the
        chromosome, typically promoters
    motif_list : Iterable[str]
        An Iterable containing the names of transcription factors of
        interest
    chrom : str
        The name of the chromosome of interest
    start_ind : int
        The starting index of the region DataFrame (bed_df)
    final_ind : int
        The final index of the region DataFrame (bed_df)
    score_threshold : float, optional
        The score required for selection, by default 400

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered edges from the 
        regions of interest
    """
    
    df = pd.DataFrame()
    
    for transcript in bed_df.index[start_ind:final_ind]:
        # Retrieve the start and end points
        start,end = bed_df.loc[transcript][['start', 'end']]
        # Load all matches in that region from the bigbed file
        motifs = bioframe.read_bigbed(bb_ref, chrom, start=start, end=end)
        # Filter only the transcription factors in the list
        max_scores = motifs[motifs['TFName'].isin(motif_list)].groupby(
            'TFName')[['score', 'name']].max()
        # Filter only high enough scores
        max_scores = max_scores[max_scores['score'] >= score_threshold]
        max_scores.reset_index(inplace=True)
        # Add the transcript (region) name for identification
        max_scores['transcript'] = transcript
        # Append the results to the dataframe
        df = pd.concat((df, max_scores), ignore_index=True, copy=False)

    return df


def filter_edges_helper(
    input_tuple: FILTER_INPUT
) -> pd.DataFrame:
    """
    Serves as a wrapper around the filter_edges function that only
    takes a single input, allowing it to be processed by map_async from
    the multiprocessing module.

    Parameters
    ----------
    input_tuple : FILTER_INPUT
        A tuple of all inputs to the filter_edges function, for more
        details refer to its docstring

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the filtered edges corresponding
        to the given input
    """

    return filter_edges(*input_tuple)