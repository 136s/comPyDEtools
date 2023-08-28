#!/usr/bin/env python

import pickle

from . import parser
from .condition import set_condition, CONDITION
from .core import Paper
from .utils import run_commands


def main():
    # generate datasets
    paper = Paper(nrep=CONDITION.nrep)
    paper.generate_datasets()
    # run analysis
    for anal_res in run_commands(CONDITION.analysis.cmds):
        print(anal_res)
    # calc metrics
    paper.calc_metrics()
    paper.res2d.to_csv("paper.csv")
    # draw figures
    paper.draw()
    with open("paper.pickle", "wb") as f:
        pickle.dump(paper, f)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.condition:
        for cnd in args.condition:
            set_condition(cnd)
        main()
