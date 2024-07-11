### Imports ###
import bioframe

import pandas as pd

from pathlib import Path
from typing import Iterable, Tuple

FILTER_INPUT = Tuple[str, pd.DataFrame, Iterable[str], str, int, int, float]

### Functions ###
def filter_edges(
    bb_ref: Path, 
    bed_df: pd.DataFrame, 
    motif_list: Iterable[str], 
    chrom: str, 
    start_ind: int, 
    final_ind: int, 
    score_threshold: float = 400,
) -> pd.DataFrame:
    """
    Filters possible binding site matches for the provided transcription
    factors from a bigbed file. This is done for a number of continuous
    regions of a single chromosome subject to a score threshold.

    Parameters
    ----------
    bb_ref : Path
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