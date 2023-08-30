#!/usr/bin/env python

import pickle
import sys

from . import parser
from .condition import set_condition, CONDITION, COMP_RES_DIR
from .core import Paper
from .utils import run_commands


def main(just_draw: bool = False):
    paper = Paper(nrep=CONDITION.nrep)
    if not just_draw:
        # generate datasets
        print("Generating datasets...", file=sys.stdout)
        paper.generate_datasets()
        # run analysis
        print("Running DE analysis...", file=sys.stdout)
        for anal_res in run_commands(CONDITION.analysis.cmds):
            print(anal_res, file=sys.stdout)
        # calc metrics
        print("Calculating metrics of DE analysis...", file=sys.stdout)
        paper.calc_metrics()
        # save metrics values
        COMP_RES_DIR.mkdir(exist_ok=True, parents=True)
        paper.res2d.to_csv(
            COMP_RES_DIR.joinpath("metrics_values.csv"),
            lineterminator="\n",
            encoding="utf-8-sig",
        )
    # draw figures
    print("Drawing figures...", file=sys.stdout)
    paper.make()
    with open(COMP_RES_DIR.joinpath("paper.pickle"), "wb") as f:
        pickle.dump(paper, f)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.condition:
        for cnd in args.condition:
            set_condition(cnd)
        main(just_draw=args.draw)
