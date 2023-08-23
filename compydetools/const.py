#!/usr/bin/env python

"""
Constant classes
"""

from enum import Enum as OriginalEnum, unique
from pathlib import Path
from typing import Any, Self


# 文字型と Path 型の Union 型
StrPath = str | Path


class Enum(OriginalEnum):
    @classmethod
    @property
    def values(cls) -> list[Any]:
        return [i.value for i in cls]


@unique
class Simul(Enum):
    KIRC = "KIRC"
    Bottomly = "Bottomly"
    mKdB = "mKdB"
    mBdK = "mBdK"

    @classmethod
    def ngenes(cls, simul_data: Self | str) -> int:
        ngenes_dict: dict[Self, int] = {
            cls.KIRC.value: 10000,
            cls.Bottomly.value: 5000,
            cls.mKdB.value: 5000,
            cls.mBdK.value: 5000,
        }
        if not isinstance(simul_data, str):
            simul_data = simul_data.value
        return ngenes_dict.get(simul_data)


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
    SEED: int = 368697996
