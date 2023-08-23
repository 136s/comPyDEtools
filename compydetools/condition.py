#!/usr/bin/env python

"""
Handle simulation conditions.
"""

from pathlib import Path
from pprint import pformat

from hydra import compose, initialize
from omegaconf import DictConfig

from .const import Simul, Disp, Outlier, Metrics, Method, Default


CONDITION = DictConfig(
    {
        "simul_data": Simul.values,
        "disp_type": Disp.values,
        "frac_up": Default.FRAC_UP,
        "nsample": Default.NSAMPLE,
        "outlier_mode": Outlier.values,
        "pde": Default.PDE,
        "metrics_type": Metrics.values,
        "method_type": Method.values,
        "nrep": Default.NREP,
    }
)


def set_condition(condition_path: str | Path) -> None:
    condition_path = Path(condition_path)
    if condition_path.is_file():
        global CONDITION
        with initialize(
            version_base=None, config_path=condition_path.parent.as_posix()
        ):
            overriding = compose(config_name=condition_path.stem)
        CONDITION.merge_with(overriding)
    print(f"CINDITION: {pformat(dict(CONDITION))}")
