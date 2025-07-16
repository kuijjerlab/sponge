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
import os
import shutil

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

from .sponge import Sponge

### Functions ###
def cli(
) -> None:
    """
    Command line interface to the SPONGE class. Execute with the --help
    option for more details.
    """

    DESCRIPTION = """
    SPONGE - Simple Prior Omics Network GEnerator.
    Generates prior motif and PPI networks, usable by other NetZoo tools
    (most notably PANDA).
    Uses the Ensembl, JASPAR, UniProt, NCBI, and STRING databases.
    Developed by Ladislav Hovan (ladislav.hovan@ncmbm.uio.no).
    """
    EPILOG = """
    Code available under GPL-3.0 license:
    https://github.com/kuijjerlab/sponge
    """

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
        description=DESCRIPTION, epilog=EPILOG)

    parser.add_argument('-t', '--temp-folder', dest='temp_folder',
        help='folder to save temporary files to',
        default='.sponge_temp', metavar='DIR')
    parser.add_argument('-c', '--config', dest='config_file',
        help='file containing the configuration',
        default='user_config.yaml', metavar='FILE')
    parser.add_argument('-e', '--example', dest='example_config',
        help='create an example config file with the default values called '
            'user_config.yaml and exit',
        action='store_true')

    args = parser.parse_args()

    if args.example_config:
        file_dir = Path(__file__).parents[0]
        shutil.copy(os.path.join(file_dir, 'user_config.yaml'),
            'user_config.yaml')
    else:
        _ = Sponge(
            temp_folder=args.temp_folder,
            config=args.config_file,
        )