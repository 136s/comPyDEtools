#!/usr/bin/env python

"""
Handle simulation conditions.
"""

from pathlib import Path
from pprint import pformat

from hydra import compose, initialize_config_dir
from omegaconf import DictConfig, OmegaConf

from .const import Enum, StrPath, Simul, Disp, Outlier, Metrics, Method, Default


CONDTYPE_CALSS: dict[str, Enum] = {
    "simul_data": Simul,
    "disp_type": Disp,
    "outlier_mode": Outlier,
    "metrics_type": Metrics,
    "method_type": Method,
}

CONDITION = DictConfig(
    {cond_type: list(enum_class) for cond_type, enum_class in CONDTYPE_CALSS.items()}
    | {
        "frac_up": Default.FRAC_UP,
        "nsample": Default.NSAMPLE,
        "pde": Default.PDE,
        "nrep": Default.NREP,
        "dirs": {
            "de_input": "input",
            "de_output": "output",
            "compre_result": "result",
        },
    }
)

DE_INPUT_DIR = Path(CONDITION.dirs.de_input)
DE_OUTPUT_DIR = Path(CONDITION.dirs.de_output)
COMP_RES_DIR = Path(CONDITION.dirs.compre_result)


def set_condition(condition_path: StrPath) -> None:
    condition_path = Path(condition_path).resolve()
    if condition_path.is_file():
        global CONDITION
        with initialize_config_dir(
            version_base=None, config_dir=condition_path.parent.as_posix()
        ):
            overriding = compose(config_name=condition_path.stem)
            # cast from str to Enum.
            OmegaConf.set_struct(overriding, False)
            for cond_type, enum_class in CONDTYPE_CALSS.items():
                overriding[cond_type] = [enum_class[i] for i in overriding[cond_type]]
            OmegaConf.set_struct(overriding, True)
        CONDITION.merge_with(overriding)
    print(f"CINDITION: {pformat(dict(CONDITION))}")
    global DE_INPUT_DIR
    DE_INPUT_DIR = Path(CONDITION.dirs.de_input)
    global DE_OUTPUT_DIR
    DE_OUTPUT_DIR = Path(CONDITION.dirs.de_output)
    global COMP_RES_DIR
    COMP_RES_DIR = Path(CONDITION.dirs.compre_result)
