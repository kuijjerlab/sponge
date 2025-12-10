# Copyright (C) 2025 Ladislav Hovan <ladislav.hovan@ncmbm.uio.no>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public License as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public License along
# with this library. If not, see <https://www.gnu.org/licenses/>.

### Imports ###
import gzip
import requests

import pandas as pd
import xml.etree.ElementTree as et

from io import BytesIO
from pathlib import Path
from tqdm import tqdm
from typing import Iterable, List, Optional, Tuple, Union

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


def create_xml_query(
    dataset_name: str,
    requested_fields: Iterable[str],
) -> str:
    """
    Formulates an XML query to retrieve specified field from a dataset
    using Ensembl BioMart.

    Parameters
    ----------
    dataset_name : str
        Name of the dataset
    requested_fields : Iterable[str]
        Fields to be retrieved

    Returns
    -------
    str
        Formulated XML query
    """

    # Build up the XML query
    xml_query = et.Element('Query', attrib=dict(virtualSchemaName='default',
        formatter='TSV', header='1', uniqueRows='0', count='',
        datasetConfigVersion='0.6'))
    dataset = et.SubElement(xml_query, 'Dataset',
        attrib=dict(name=dataset_name, interface='default'))
    for field in requested_fields:
        _ = et.SubElement(dataset, 'Attribute', attrib=dict(name=field))
    # Convert to a string with a declaration
    query_string = et.tostring(xml_query, xml_declaration=True,
        encoding='unicode')

    return query_string


def retrieve_ensembl_data(
    dataset_name: str,
    requested_fields: Iterable[str],
    ensembl_url: str,
) -> BytesIO:
    """
    Retrieves specified fields from an Ensembl dataset by querying
    BioMart.

    Parameters
    ----------
    dataset_name : str
        Name of the dataset
    requested_fields : Iterable[str]
        Fields to be retrieved
    ensembl_url : str
        URL for BioMart

    Returns
    -------
    BytesIO
        Bytes retrieved from the server
    """

    xml_query = create_xml_query(dataset_name, requested_fields)
    REQUEST_STRING = '/martservice?query='
    link = ensembl_url + REQUEST_STRING + xml_query
    r = requests.get(link, stream=True)
    r.raise_for_status()
    bytes = download_with_progress(r)

    return bytes


def get_ensembl_version(
    ensembl_rest: str,
) -> str:
    """
    Returns the full version of the genome assembly used by the
    Ensembl server (e.g. GRCh38).

    Parameters
    ----------
    ensembl_rest : str
        URL for the REST interface of Ensembl

    Returns
    -------
    str
        Full version of the genome assembly used by Ensembl
    """

    # Request the assembly information from Ensembl REST
    REQUEST_STRING = "/info/assembly/homo_sapiens?"
    r = requests.get(ensembl_rest + REQUEST_STRING,
        headers={ "Content-Type" : "application/json"})
    r.raise_for_status()
    decoded = r.json()

    return decoded['assembly_name']


def get_chromosome_mapping(
    assembly: str,
    mapping_url: str,
) -> pd.Series:
    """
    Returns a  pandas Series which can be used to map Ensembl chromosome
    names to UCSC for a provided genome assembly. If it is not
    recognised, None is returned.

    Parameters
    ----------
    assembly : str
        Assembly for which to provide the mapping (e.g. hg38)
    mapping_url : str
        URL to retrieve chromosome mapping from

    Returns
    -------
    Tuple[pd.Series, pd.Series]
        Pandas Series providing chromosome name mapping from Ensembl
        to UCSC
    """

    if assembly[:2] == 'hg':
        # The mapping can maybe be retrieved from a chromAlias.tsv file
        print (f'Retrieving chromosome name mapping for {assembly}...')
        try:
            f = gzip.open(download_with_progress(mapping_url.format(
                genome_assembly=assembly)))
            header_fields = ['alt', 'ucsc', 'notes']
            chrom_df = pd.read_csv(f, sep='\t', names=header_fields)
            # This mapping is unambiguous even if things other than Ensembl
            # are included in the index
            return chrom_df.set_index('alt')['ucsc']
        except requests.exceptions.HTTPError:
            print (f'Failed to retrieve mapping for the assembly {assembly}.')
            return None
    else:
        print ('No chromosome name mapping available for the assembly',
            f'{assembly}.')
        return None