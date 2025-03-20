### Imports ###
import pandas as pd

from pathlib import Path
from typing import Iterable

### Class definition ###
class FileWriter:
    # Functions
    def __init__(
        self,
    ):
        
        pass


    def write_network_file(
        self,
        df: pd.DataFrame,
        node_columns: Iterable[str],
        weight_column: str,
        file_name: Path,
    ):
        
        sorted_df = df.sort_values(by=node_columns)
        sorted_df[node_columns + [weight_column]].to_csv(
            file_name, sep='\t', index=False, header=False)