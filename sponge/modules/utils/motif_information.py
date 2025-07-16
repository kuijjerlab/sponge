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
import pandas as pd

from Bio.motifs.jaspar import Motif
from math import log2

### Functions ###
def plogp(
    x: float,
) -> float:
    """
    Returns x*log2(x) for a number, handles the 0 case properly.

    Parameters
    ----------
    x : float
        Input value

    Returns
    -------
    float
        Value of x*log2(x)
    """

    if x == 0:
        return 0
    else:
        return x*log2(x)


def calculate_ic(
    motif: Motif,
) -> float:
    """
    Calculates the information content for a given motif, assuming equal
    ACGT distribution.

    Parameters
    ----------
    motif : Motif
        JASPAR Motif object

    Returns
    -------
    float
        Information content of the motif
    """

    df = pd.DataFrame(motif.pwm)
    # Calculate the IC for each position in the motif
    df['IC'] = df.apply(lambda x: 2 + sum([plogp(x[y]) for y in
        ['A', 'C', 'G', 'T']]), axis=1)

    # Return the total IC for the whole motif
    return df['IC'].sum()