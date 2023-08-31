#!/usr/bin/env python

from pathlib import Path
import pickle
import sys

from . import parser
from .condition import set_condition, CONDITION
from .core import Paper
from .utils import run_commands


def main(just_generate: bool = False, just_draw: bool = False):
    comp_res_dir = Path(CONDITION.dirs.result)
    concat_condition = "_".join(
        [
            "{" + "|".join([i.name for i in CONDITION.simul_data]) + "}",
            "{" + "|".join([i.name for i in CONDITION.disp_type]) + "}",
            "upFrac{" + "|".join([str(i) for i in CONDITION.frac_up]) + "}",
            "{" + "|".join([str(i) for i in CONDITION.nsample]) + "}spc",
            "{" + "|".join([i.name for i in CONDITION.outlier_mode]) + "}",
            "{" + "|".join([str(i) for i in CONDITION.pde]) + "}DE",
            f"nrep{CONDITION.nrep}",
        ]
    )
    paper = Paper(nrep=CONDITION.nrep)
    if not just_draw:
        # generate datasets
        print("Generating datasets...", file=sys.stdout)
        paper.generate_datasets()
    if not (just_draw or just_generate):
        # run analysis
        print("Running DE analysis...", file=sys.stdout)
        for anal_res in run_commands(CONDITION.analysis.cmds):
            print(anal_res, file=sys.stdout)
        # calc metrics
        print("Calculating metrics of DE analysis...", file=sys.stdout)
        paper.calc_metrics()
        # save metrics values
        comp_res_dir.mkdir(exist_ok=True, parents=True)
        paper.res2d.to_csv(
            comp_res_dir.joinpath(concat_condition + ".csv"),
            lineterminator="\n",
            encoding="utf-8-sig",
        )
    if not just_generate:
        # draw figures
        print("Drawing figures...", file=sys.stdout)
        paper.make()
        with open(comp_res_dir.joinpath(concat_condition + ".pickle"), "wb") as f:
            pickle.dump(paper, f)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.condition:
        for cnd in args.condition:
            set_condition(cnd)
        main(just_generate=args.generate, just_draw=args.draw)
