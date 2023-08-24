#!/usr/bin/env python

"""
Constant classes
"""

from enum import Enum as OriginalEnum, unique
from pathlib import Path
from typing import Self


StrPath = str | Path
DE_INPUT_DIR = Path("input")
DE_OUTPUT_DIR = Path("output")
COMP_RES_DIR = Path("result")


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
    cutoff = "Best threshold"
    f1score = "F1-score"
    kappa = "Cohen's kappa"


@unique
class Method(Enum):
    fc = "Fold change"
    nc = "Normalized change"
    rp = "Rank product"
    cp = "Combined probability"
    deseq2 = "DEseq2"


class Default:
    FRAC_UP: tuple[float] = (0.5, 0.7, 0.9)
    NSAMPLE: tuple[int] = (3, 10)
    PDE: tuple[float] = (0.27, 5, 10, 30, 60)
    NREP: int = 50
    SEED: int = 330252033
