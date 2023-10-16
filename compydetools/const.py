#!/usr/bin/env python

"""
Constant classes
"""

from enum import Enum as OriginalEnum, unique
from pathlib import Path
from typing import Literal, Self
import warnings


StrPath = str | Path
MetricsInput = dict[Literal["y_true", "y_score", "y_pred"], list[float]]


class Enum(OriginalEnum):
    def __repr__(self) -> str:
        return str(self)

    @classmethod
    @property
    def first(cls):
        return list(cls)[0]


@unique
class Simul(Enum):
    KIRC = "KIRC"
    Bottomly = "Bottomly"
    mKdB = "mKdB"
    mBdK = "mBdK"

    @classmethod
    def ngenes(cls, simul_data: Self) -> int:
        match simul_data:
            case cls.KIRC:
                return 10000
            case cls.Bottomly | cls.mKdB | cls.mBdK:
                return 5000
            case _:
                raise AttributeError


@unique
class Disp(Enum):
    same = "same"
    different = "different"


@unique
class Outlier(Enum):
    D = "no random outlier"
    R = "random outlier test"
    OS = "outlying dispersion sample contained"
    DL = "lowered dispersion test"


@unique
class Metrics(Enum):
    auc = "AUC"
    tpr = "TPR"
    fdr = "True FDR"
    f1score = "F1-score"
    kappa = "Cohen's kappa"
    cutoff = "Best threshold"


@unique
class Method(Enum):
    deseq2 = "DEseq2"


def update_method_class(d: dict):
    try:
        global Method
        Method = Enum("Method", d)
    except ValueError:
        warnings.warn(f"Invalid methods = {d}")


class Default:
    FRAC_UP: tuple[float] = (0.5, 0.7, 0.9)
    NSAMPLE: tuple[int] = (3, 10)
    PDE: tuple[float] = (0.27, 5, 10, 30, 60)
    NREP: int = 50
    SEED: int = 330252033
    PALETTE: list[str] = [
        "#ff4b00",
        "#fff100",
        "#03af7a",
        "#ed7ca2",
        "#005aff",
        "#f6aa00",
        "#990099",
    ]
