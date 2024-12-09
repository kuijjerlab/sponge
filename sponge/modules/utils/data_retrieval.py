### Imports ###
import requests

from io import BytesIO
from pathlib import Path
from tqdm import tqdm
from typing import List, Optional, Union

### Functions ###
def download_with_progress(
    url: Union[List[str], str, requests.models.Response],
    file_path: Optional[Path] = None,
    desc: str = 'response',
) -> Optional[BytesIO]:
    """
    Downloads from a given URL or retrieves a response to a given
    request while providing a progress bar.

    Parameters
    ----------
    url : Union[str, requests.models.Response]
        URL or response to be processed
    file_path : Optional[Path], optional
        File path for saving or None to save into a BytesIO object,
        by default None
    desc : str, optional
        Description to show, by default 'response'

    Returns
    -------
    Optional[BytesIO]
        BytesIO object containing the data or None if file_path was
        not set to None
    """

    # Determine the type of request
    if type(url) == str:
        try:
            request = requests.get(url, stream=True)
        except requests.exceptions.SSLError as ssl:
            print ('The following verification error has occured:')
            print (ssl)
            print ('Retrying without verification')
            request = requests.get(url, stream=True, verify=False)
        # Client or server errors
        request.raise_for_status()
    elif isinstance(url, List):
        # Multiple possible URLs, use the first one that works
        for pos,u in enumerate(url):
            try:
                return download_with_progress(u,
                    file_path=file_path, desc=desc)
            except requests.exceptions.ConnectionError as conn:
                if pos < len(url) - 1:
                    print ('The following URL was unreachable:')
                    print (u)
                    print ('Trying the next one')
                else:
                    raise conn
            except requests.exceptions.HTTPError as http:
                if pos < len(url) - 1:
                    print ('An HTTP error was raised when connecting to this '
                        'URL:')
                    print (u)
                    print ('Trying the next one')
                else:
                    raise http
    else:
        request = url
    total = int(request.headers.get('content-length', 0))
    # Determine whether to save data to a file or object
    if file_path is None:
        stream = BytesIO()
    else:
        stream = open(file_path, 'wb')
        desc = file_path

    # Download with a progress bar using tqdm
    with tqdm(desc=desc, total=total, unit='iB', unit_scale=True,
        unit_divisor=1024) as bar:
        for data in request.iter_content(chunk_size=1024):
            size = stream.write(data)
            bar.update(size)

    if file_path is None:
        return BytesIO(stream.getvalue())