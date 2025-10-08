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