### Imports ###
from pyjaspar import jaspar_releases, JASPAR_LATEST_RELEASE
from typing import Optional

### Functions ###
def process_jaspar_version(
    jaspar_release: Optional[str],
) -> str:
    """
    Processes a desired JASPAR release by attempting to match it to
    an existing one. If None is provided, will return the newest
    JASPAR release.

    Parameters
    ----------
    jaspar_release : Optional[str]
        Provided desired JASPAR release or None to get the newest

    Returns
    -------
    str
        Processed JASPAR release that matches the provided one

    Raises
    ------
    ValueError
        If the provided JASPAR release couldn't be matched to an
        existing one
    """

    if jaspar_release is None:
        # Just keep the current object, log the release version
        jaspar_release = JASPAR_LATEST_RELEASE
    else:
        jaspar_available = [jr for jr in jaspar_releases.keys()]
        alt_name = 'JASPAR' + jaspar_release
        if jaspar_release in jaspar_available:
            # Release found as specified
            pass
        elif alt_name in jaspar_available:
            # Try adding JASPAR to the provided release
            # Converts e.g. 2022 to JASPAR2022 (actual release name)
            print (f'Found {alt_name} in available releases, assuming '
                'this matches your choice')
            jaspar_release = alt_name
        else:
            error_str = ('The specified JASPAR release '
                f'({jaspar_release}) is not available.\n'
                'Available releases: ' +
                ', '.join(jaspar_available))
            raise ValueError(error_str)

    return jaspar_release