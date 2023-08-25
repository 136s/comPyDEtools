#!/usr/bin/env python

import subprocess

from . import parser
from .condition import set_condition, CONDITION
from .core import AllPages


def main():
    # generate datasets
    allpages = AllPages()
    allpages.generate_datasets()
    # run analysis
    for cmd in CONDITION.analysis.cmds:
        subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.condition:
        for cnd in args.condition:
            set_condition(cnd)
        main()
