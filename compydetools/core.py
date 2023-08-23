#!/usr/bin/env python

"""
Simulation core classes
"""

from dataclasses import dataclass, field

from .condition import CONDITION
from .const import Simul, Disp, Outlier, Metrics, Method, Default


@dataclass
class Result:
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.values[0])
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.values[0])
    method_type: Method = field(default=Method.values[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    value: float = field(init=False, repr=False)


@dataclass
class Dataset:
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.values[0])
    pde: float = field(default=Default.PDE[0])
    seed: int = field(default=Default.SEED)


@dataclass
class Box:
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.values[0])
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.values[0])
    method_type: Method = field(default=Method.values[0])
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
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.values[0])
    pde: float = field(default=Default.PDE[0])
    metrics_type: Metrics = field(default=Metrics.values[0])
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
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
    frac_up: float = field(default=Default.FRAC_UP[0])
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.values[0])
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
    simul_data: Simul = field(default=Simul.values[0])
    disp_type: Disp = field(default=Disp.values[0])
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
