#!/usr/bin/env python

"""
Simulation core classes
"""

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .condition import CONDITION
from .const import (
    DE_INPUT_DIR,
    Simul,
    Disp,
    Outlier,
    Metrics,
    Method,
    Default,
)
from .generation import synthetic_data_simulation


@dataclass
class Result:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.first)
    method_type: Method = field(default=Method.first)
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    value: float = field(init=False, repr=False)


@dataclass
class Dataset:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    seed: int = field(default=Default.SEED)
    ngenes: int = field(init=False, repr=False)
    nde: int = field(init=False, repr=False)
    counts: pd.DataFrame = field(init=False, repr=False)
    cond_str: str = field(init=False, repr=False)
    countpath: Path = field(init=False, repr=False)
    configpath: Path = field(init=False, repr=False)

    def convert_pde(self) -> None:
        self.ngenes = Simul.ngenes(self.simul_data)
        self.nde = round(self.ngenes * self.pde / 100)

    def generate(self) -> pd.DataFrame:
        if not hasattr(self, "counts"):
            self.convert_pde()
            self.counts = synthetic_data_simulation(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                ngenes=self.ngenes,
                nde=self.nde,
                seed=self.seed,
            )
        return self.counts

    def set_path(self):
        self.cond_str = "_".join(
            [
                self.simul_data.name,
                self.disp_type.name,
                f"upFrac{self.frac_up}",
                f"{self.nsample}spc",
                self.outlier_mode.name,
                f"{self.nde}DE",
            ]
        )
        self.countpath = DE_INPUT_DIR.joinpath(
            self.cond_str, f"{self.cond_str}_rep{self.seed}.tsv"
        )
        self.configpath = self.countpath.with_suffix(".yaml")

    def save(self):
        self.generate()
        self.set_path()
        self.countpath.parent.mkdir(exist_ok=True, parents=True)
        self.counts.to_csv(
            self.countpath, sep="\t", lineterminator="\n", encoding="utf-8"
        )


@dataclass
class Box:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.first)
    method_type: Method = field(default=Method.first)
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    datasets: list[Dataset] = field(default_factory=list, init=False, repr=False)
    results: list[Result] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        for seed_ in range(self.seed, self.seed + self.nrep):
            dataset = Dataset(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                pde=self.pde,
                seed=seed_,
            )
            self.datasets.append(dataset)


@dataclass
class Plot:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.first)
    boxes: list[Box] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for method_type in CONDITION.method_type:
            box = Box(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                pde=self.pde,
                metrics_type=self.metrics_type,
                method_type=method_type,
                nrep=CONDITION.nrep,
            )
            self.boxes.append(box)


@dataclass
class Figure:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    plots: list[Plot] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for pde in CONDITION.pde:
            for metrics_type in CONDITION.metrics_type:
                plot = Plot(
                    simul_data=self.simul_data,
                    disp_type=self.disp_type,
                    frac_up=self.frac_up,
                    nsample=self.nsample,
                    outlier_mode=self.outlier_mode,
                    pde=pde,
                    metrics_type=metrics_type,
                )
                self.plots.append(plot)


@dataclass
class Page:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    figures: list[Figure] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for nsample in CONDITION.nsample:
            for outlier_mode in CONDITION.outlier_mode:
                figure = Figure(
                    simul_data=self.simul_data,
                    disp_type=self.disp_type,
                    frac_up=self.frac_up,
                    nsample=nsample,
                    outlier_mode=outlier_mode,
                )
                self.figures.append(figure)
