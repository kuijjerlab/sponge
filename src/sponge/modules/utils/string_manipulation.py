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