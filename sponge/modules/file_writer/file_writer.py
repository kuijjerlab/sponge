### Imports ###
import pandas as pd

from pathlib import Path
from typing import Iterable

### Class definition ###
class FileWriter:
    # Methods
    def __init__(
        self,
    ):
        """
        Class which writes prior networks to files.
        """
        
        pass


    def write_network_file(
        self,
        df: pd.DataFrame,
        node_columns: Iterable[str],
        weight_column: str,
        file_name: Path,
    ):
        """
        Writes the selected columns of a DataFrame to the provided path.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing the network
        node_columns : Iterable[str]
            Columns to be used as annotations for the edges
        weight_column : str
            Column with weights for the network
        file_name : Path
            Where to write the network to
        """
        
        sorted_df = df.sort_values(by=node_columns)
        sorted_df[node_columns + [weight_column]].to_csv(
            file_name, sep='\t', index=False, header=False)