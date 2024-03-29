#!/usr/bin/env python

"""
Handle simulation conditions.
"""

from pathlib import Path
from pprint import pformat

from hydra import compose, initialize_config_dir
from omegaconf import DictConfig, ListConfig, OmegaConf

from .const import (
    Enum,
    StrPath,
    Simul,
    Disp,
    Outlier,
    Metrics,
    Method,
    Default,
    update_method_class,
)
from .generation import GENE_DESCRIPTION_COLNAME


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
        "analysis": {
            "cmds": [],
            "res": "",
            "de_true": GENE_DESCRIPTION_COLNAME,
            "de_score": "padj",
            "de_score_threshold": 0.1,
        },
        "dirs": {
            "dataset": "input",
            "result": "result",
        },
    }
)


def set_condition(condition_path: StrPath) -> None:
    condition_path = Path(condition_path).resolve()
    if condition_path.is_file():
        global CONDITION
        with initialize_config_dir(
            version_base=None, config_dir=condition_path.parent.as_posix()
        ):
            overriding = compose(config_name=condition_path.stem)
            # apply user methods to `const.Method`
            if isinstance(
                user_method_types := overriding.get("method_type"), DictConfig
            ):
                update_method_class(user_method_types)
                # from dict to list
                overriding.method_type = ListConfig(list(overriding.method_type.keys()))
            # cast from str to Enum.
            OmegaConf.set_struct(overriding, False)
            for cond_type, enum_cls in CONDTYPE_CALSS.items():
                if cond_type == "method_type":
                    from .const import Method  # reload Method after updating

                    overriding[cond_type] = [Method[i] for i in overriding[cond_type]]
                else:
                    overriding[cond_type] = [enum_cls[i] for i in overriding[cond_type]]
            OmegaConf.set_struct(overriding, True)
        CONDITION.merge_with(overriding)
    print(f"CINDITION: {pformat(dict(CONDITION), sort_dicts=False)}")
