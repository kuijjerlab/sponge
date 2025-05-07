### Imports ###
from typing import Mapping

### Functions ###
def recursive_update(
    dictionary: dict,
    update: dict,
) -> dict:
    """
    Recursively updates the provided dictionary with the
    provided update. This makes multiple level updates possible,
    as the changes will be added, but they will not overwrite.

    Parameters
    ----------
    dictionary : dict
        Dictionary to update
    update : dict
        Dictionary to update with

    Returns
    -------
    dict
        Updated dictionary
    """

    for k,v in update.items():
        if isinstance(dictionary, Mapping):
            if isinstance(v, Mapping):
                r = recursive_update(dictionary.get(k, {}), v)
                dictionary[k] = r
            else:
                dictionary[k] = update[k]
        else:
            dictionary = {k: update[k]}

    return dictionary