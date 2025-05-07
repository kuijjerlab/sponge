### Imports ###
from datetime import datetime
from typing import Union

### Functions ###
def adjust_gene_name(
    gene: str,
) -> str:
    """
    Adjusts the gene name by converting the last two letters to
    lowercase. This is typically done to enhance name matching. Deals
    with dimers properly.

    Parameters
    ----------
    gene : str
        Provided gene name

    Returns
    -------
    str
        Adjusted gene name
    """

    genes = gene.split('::')

    return '::'.join([g[:-2] + g[-2:].lower() for g in genes])


def parse_datetime(
    datetime: Union[str, datetime],
) -> str:
    """
    Converts the provided datetime object into a formatted string,
    or returns the provided string.

    Parameters
    ----------
    datetime : Union[str, datetime]
        Provided string or datetime object

    Returns
    -------
    str
        Provided string or provided datetime expressed as a string
    """

    if type(datetime) == str:
        return datetime
    else:
        return datetime.strftime('%d/%m/%Y, %H:%M:%S')