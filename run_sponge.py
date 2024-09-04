#!/usr/bin/env python

import argparse

from sponge import Sponge

if __name__ == '__main__':
    sponge_obj = Sponge(
        run_default=True,

    )
    sponge_obj.show_fingerprint()