#!/usr/bin/env python

"""
Simulation core classes
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
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
from .generation import synthetic_data_simulation, GENE_ID_COLNAME


@dataclass
class Result:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    seed: int = field(default=Default.SEED)
    metrics_type: Metrics = field(default=Metrics.first)
    method_type: Method = field(default=Method.first)
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
    filepath: Path = field(init=False, repr=False)

    def convert_pde(self) -> None:
        self.ngenes = Simul.ngenes(self.simul_data)
        self.nde = round(self.ngenes * self.pde / 100)

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
        self.filepath = DE_INPUT_DIR.joinpath(
            self.cond_str, f"{self.cond_str}_rep{self.seed}.tsv"
        )

    def generate(self) -> pd.DataFrame:
        if not hasattr(self, "counts"):
            self.convert_pde()
            self.set_path()
            if self.filepath.is_file():
                self.counts = pd.read_table(self.filepath, index_col=GENE_ID_COLNAME)
            else:
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
                self.filepath.parent.mkdir(exist_ok=True, parents=True)
                self.counts.to_csv(
                    self.filepath, sep="\t", lineterminator="\n", encoding="utf-8"
                )
        return self.counts


@dataclass
class Box:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.first)
    method_type: Method = field(default=Method.first)
    results: list[Result] = field(init=False, repr=False)


@dataclass
class Plot:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.first)
    boxes: list[Box] = field(default_factory=list, init=False, repr=False)
    datasets: list[Dataset] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for method_type in CONDITION.method_type:
            box = Box(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nrep=self.nrep,
                seed=self.seed,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                pde=self.pde,
                metrics_type=self.metrics_type,
                method_type=method_type,
            )
            self.boxes.append(box)
        ds_seeds = np.random.default_rng(self.seed).integers(
            0, np.iinfo(np.int64).max, size=self.nrep
        )
        for ds_seed in ds_seeds:
            dataset = Dataset(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                pde=self.pde,
                seed=ds_seed,
            )
            self.datasets.append(dataset)


@dataclass
class Figure:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
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
                    nrep=self.nrep,
                    seed=self.seed,
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
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    figures: list[Figure] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for nsample in CONDITION.nsample:
            for outlier_mode in CONDITION.outlier_mode:
                figure = Figure(
                    simul_data=self.simul_data,
                    disp_type=self.disp_type,
                    frac_up=self.frac_up,
                    nrep=self.nrep,
                    seed=self.seed,
                    nsample=nsample,
                    outlier_mode=outlier_mode,
                )
                self.figures.append(figure)
