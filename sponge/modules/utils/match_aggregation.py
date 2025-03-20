### Imports ###
import numpy as np
import pandas as pd

from pathlib import Path

### Functions ###
def aggregate_matches(
    match_df: pd.DataFrame,
    ensembl_file: Path,
    animal_to_human: dict,
    use_gene_names: bool,
    protein_coding_only: bool,
) -> pd.DataFrame:
    """
    Aggregates all the matches corresponding to individual
    transcripts into genes, creating a transcription factor - gene
    matrix. Stores the result internally.

    Parameters
    ----------
    ensembl_file : Optional[Path], optional
        Path to an Ensembl file or None to use cache or
        download it, by default None
    prompt : bool, optional
        Whether to prompt before downloading, by default True
    use_gene_names : Optional[bool], optional
        Whether to use gene names instead of Ensembl IDs or None
        to follow the option from the initialisation,
        by default None
    protein_coding_only : Optional[bool], optional
        Whether to restrict the selection to only protein coding
        genes or None to follow the option from the initialisation,
        by default None
    """

    ensembl = pd.read_csv(ensembl_file, sep='\t')

    # Add the Ensembl data (gene names) to the edges previously found
    motif_df = match_df.join(other=ensembl.set_index(
        'Transcript stable ID'), on='transcript')
    print ()
    print ('Number of TF - transcript edges:', len(motif_df))
    if protein_coding_only:
        motif_df = motif_df[motif_df['Gene type'] ==
            'protein_coding'].copy()
    # Drop columns that are not required anymore
    motif_df.drop(columns=['Gene type', 'name'], inplace=True)
    # Humanise the TF names
    motif_df['TFName'] = motif_df['TFName'].apply(lambda x:
        animal_to_human[x] if x in animal_to_human else x)
    # Ignore genes without identifiers
    motif_df.dropna(subset=['Gene stable ID'], inplace=True)
    motif_df.sort_values('score', ascending=False, inplace=True)
    # Sometimes edges are identified from multiple transcripts
    motif_df.drop_duplicates(subset=['TFName', 'Gene stable ID'],
        inplace=True)
    print ('Number of TF - gene edges:', len(motif_df))
    if use_gene_names:
        # Names are not unique - filtering needed
        # Fill empty gene names with IDs
        motif_df['Gene name'] = motif_df.apply(lambda x: x['Gene name'] if
            type(x['Gene name']) == str else x['Gene stable ID'], axis=1)
        # Count the number of edges for every name/ID pair
        name_id_matching = motif_df.groupby(
            ['Gene name', 'Gene stable ID'])['Gene name'].count()
        # Use the name for the ID that has the most edges
        id_to_name = {i[1]: i[0] for i in name_id_matching.groupby(
            level=0).idxmax().values}
        # Convert selected gene IDs to names
        motif_df['Gene name'] = motif_df['Gene stable ID'].apply(
            lambda x: id_to_name[x] if x in id_to_name else np.nan)
        # Drop the rest
        motif_df.dropna(subset='Gene name', inplace=True)
        print ('Number of TF - gene edges after name conversion:',
            len(motif_df))

    return motif_df
    