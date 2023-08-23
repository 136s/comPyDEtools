#!/usr/bin/env python

"""
Constant classes
"""

from enum import Enum as OriginalEnum, unique
from typing import Any, Self


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
    D = "D"
    R = "R"
    OS = "OS"
    DL = "DL"


@unique
class Metrics(Enum):
    auc = "auc"
    tpr = "tpr"
    fdr = "fdr"
    cutoff = "cutoff"
    f1score = "f1score"
    kappa = "kappa"


@unique
class Method(Enum):
    fc = "fc"
    nc = "nc"
    rp = "rp"
    cp = "cp"
    deseq2 = "deseq2"


class Default:
    FRAC_UP: tuple[float] = (0.5, 0.7, 0.9)
    NSAMPLE: tuple[int] = (3, 10)
    PDE: tuple[float] = (0.27, 5, 10, 30, 60)
    NREP: int = 50
    SEED: int = 368697996
