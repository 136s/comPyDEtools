#!/usr/bin/env python

"""
Simulation core classes
"""

from dataclasses import dataclass, field
from pathlib import Path

from matplotlib.figure import Figure as mplFigure
import numpy as np
import pandas as pd

from . import utils
from .condition import CONDITION, DE_INPUT_DIR, COMP_RES_DIR
from .const import Simul, Disp, Outlier, Metrics, Method, Default, MetricsInput
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
    dsid: str = field(default=None)
    method_types: list[Method] = field(default_factory=list)
    method_values: dict[Method, MetricsInput] = field(default_factory=dict, repr=False)
    metrics_types: list[Metrics] = field(default_factory=list)
    metrics_values: dict[Metrics, dict[Method, float]] = field(
        default_factory=dict, init=False, repr=False
    )
    _res2d: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.set_method_values()

    def set_method_values(self) -> None:
        if not self.method_values:
            for method_type in self.method_types:
                de_output_path = CONDITION.analysis.res.replace(
                    "{count_stem}", self.dsid
                ).replace("{method_type}", method_type.name)
                self.method_values[method_type] = utils.load_metrics_input(
                    de_output_path,
                    de_true_colname=CONDITION.analysis.de_true,
                    de_score_colname=CONDITION.analysis.de_score,
                    de_score_threshold=CONDITION.analysis.de_score_threshold,
                )

    def calc_metrics(self) -> None:
        if not self.metrics_values:
            for method_type, method_value in self.method_values.items():
                for metrics_type in self.metrics_types:
                    self.metrics_values.setdefault(metrics_type, {})
                    self.metrics_values[metrics_type][method_type] = utils.calc_metrics(
                        metrics_type, **method_value
                    )

    def set_res2d(self) -> None:
        if not hasattr(self, "_res2d"):
            self.calc_metrics()
            df = pd.DataFrame(self.metrics_values)
            df.index = df.index.map(lambda v: v.name)
            df.columns = df.columns.map(lambda v: v.name)
            self._res2d = pd.melt(
                df.reset_index(names="method"), id_vars="method", var_name="metrics"
            )
            self._res2d["simul_data"] = self.simul_data.name
            self._res2d["disp_type"] = self.disp_type.name
            self._res2d["frac_up"] = self.frac_up
            self._res2d["nsample"] = self.nsample
            self._res2d["outlier_mode"] = self.outlier_mode.name
            self._res2d["pde"] = self.pde
            self._res2d["seed"] = self.seed
            self._res2d.set_index(
                [
                    "simul_data",
                    "disp_type",
                    "frac_up",
                    "nsample",
                    "outlier_mode",
                    "pde",
                    "seed",
                    "method",
                    "metrics",
                ],
                inplace=True,
            )

    @property
    def res2d(self) -> pd.DataFrame:
        self.set_res2d()
        return self._res2d


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
    dsid: str = field(init=False, repr=False)
    filepath: Path = field(init=False, repr=False)

    def __post_init__(self) -> None:
        self.convert_pde()
        self.set_path()

    def convert_pde(self) -> None:
        self.ngenes = Simul.ngenes(self.simul_data)
        self.nde = round(self.ngenes * self.pde / 100)

    def set_path(self) -> None:
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
        self.dsid = f"{self.cond_str}_rep{self.seed}"
        self.filepath = DE_INPUT_DIR.joinpath(self.cond_str, f"{self.dsid}.tsv")

    def generate(self) -> None:
        if not hasattr(self, "counts"):
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


@dataclass
class DataPool:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    pde: float = field(default=Default.PDE[0])
    datasets: list[Dataset] = field(default_factory=list, init=False, repr=False)
    results: list[Result] = field(default_factory=list, init=False, repr=False)
    plot_values: dict[Metrics, dict[Method, list[float]]] = field(
        default_factory=dict, init=False, repr=False
    )

    def __post_init__(self) -> None:
        self.set_datasets()

    def set_datasets(self) -> None:
        ds_seeds = np.random.default_rng(self.seed).integers(
            0, np.iinfo(np.int32).max, size=self.nrep
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

    def generate_datasets(self) -> list[Dataset]:
        for dataset in self.datasets:
            dataset.generate()
        return self.datasets

    def set_results(self) -> None:
        if not self.results:
            for dataset in self.datasets:
                result = Result(
                    simul_data=dataset.simul_data,
                    disp_type=dataset.disp_type,
                    frac_up=dataset.frac_up,
                    nsample=dataset.nsample,
                    outlier_mode=dataset.outlier_mode,
                    pde=dataset.pde,
                    seed=dataset.seed,
                    dsid=dataset.dsid,
                    method_types=CONDITION.method_type,
                    metrics_types=CONDITION.metrics_type,
                )
                self.results.append(result)

    def calc_metrics(self) -> list[Result]:
        self.set_results()
        for result in self.results:
            result.calc_metrics()
        return self.results

    @property
    def res2d(self) -> pd.DataFrame:
        self.calc_metrics()
        return pd.concat([result.res2d for result in self.results])


@dataclass
class Plot:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    nsample: float = field(default=Default.NSAMPLE[0])
    outlier_mode: Outlier = field(default=Outlier.first)
    datapools: list[DataPool] = field(default_factory=list, init=False, repr=False)
    filepath: Path = field(init=False, repr=False)
    mplobj: mplFigure = field(init=False, repr=False)

    def __post_init__(self) -> None:
        for pde in CONDITION.pde:
            datapool = DataPool(
                simul_data=self.simul_data,
                disp_type=self.disp_type,
                frac_up=self.frac_up,
                nrep=self.nrep,
                seed=self.seed,
                nsample=self.nsample,
                outlier_mode=self.outlier_mode,
                pde=pde,
            )
            self.datapools.append(datapool)

    def generate_datasets(self) -> list[DataPool]:
        for datapool in self.datapools:
            datapool.generate_datasets()
        return self.datapools

    def calc_metrics(self) -> list[DataPool]:
        for datapool in self.datapools:
            datapool.calc_metrics()
        return self.datapools

    @property
    def res2d(self) -> pd.DataFrame:
        return pd.concat([result.res2d for result in self.datapools])

    def set_path(self):
        cond_str = "_".join(
            [
                self.simul_data.name,
                self.disp_type.name,
                f"upFrac{self.frac_up}",
                f"{self.nsample}spc",
                self.outlier_mode.name,
            ]
        )
        self.filepath = COMP_RES_DIR.joinpath(f"{cond_str}.svg")
        self.filepath.parent.mkdir(exist_ok=True, parents=True)

    def draw(self, save: bool = True) -> None:
        if not hasattr(self, "mplobj"):
            self.set_path()
            data = self.res2d.reset_index()
            from .const import Method  # reload Method after updating

            for colname, classname in {"method": Method, "metrics": Metrics}.items():
                data[colname] = data[colname].apply(lambda v: classname[v].value)

            plot_title = " / ".join(
                [
                    f"{self.simul_data.value} simulation",
                    f"{self.outlier_mode.value}",
                    f"{self.disp_type.value} dispersion",
                    f"SS = {self.nsample}",
                    f"Bal = {int(100*self.frac_up)}%",
                ]
            )
            self.mplobj = utils.draw_plot(
                data, plot_title=plot_title, output=self.filepath if save else None
            ).fig


@dataclass
class Figure:
    simul_data: Simul = field(default=Simul.first)
    disp_type: Disp = field(default=Disp.first)
    frac_up: float = field(default=Default.FRAC_UP[0])
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    plots: list[Plot] = field(default_factory=list, init=False, repr=False)
    filepath: Path = field(init=False, repr=False)

    def __post_init__(self) -> None:
        for nsample in CONDITION.nsample:
            for outlier_mode in CONDITION.outlier_mode:
                figure = Plot(
                    simul_data=self.simul_data,
                    disp_type=self.disp_type,
                    frac_up=self.frac_up,
                    nrep=self.nrep,
                    seed=self.seed,
                    nsample=nsample,
                    outlier_mode=outlier_mode,
                )
                self.plots.append(figure)

    def generate_datasets(self) -> list[Plot]:
        for plot in self.plots:
            plot.generate_datasets()
        return self.plots

    def calc_metrics(self) -> list[Plot]:
        for plot in self.plots:
            plot.calc_metrics()
        return self.plots

    @property
    def res2d(self) -> pd.DataFrame:
        return pd.concat([result.res2d for result in self.plots])

    def draw(self, save: bool = True) -> list[Plot]:
        for plot in self.plots:
            plot.draw(save=save)
        return self.plots

    def set_path(self):
        cond_str = "_".join(
            [self.simul_data.name, self.disp_type.name, f"upFrac{self.frac_up}"]
        )
        self.filepath = COMP_RES_DIR.joinpath(f"{cond_str}.svg")
        self.filepath.parent.mkdir(exist_ok=True, parents=True)

    def make(self, save: bool = True):
        self.draw(save=save)
        self.set_path()
        # get child plots
        children: dict[str, mplFigure] = {}
        for plot in self.plots:
            children[f"{plot.nsample}spc_{plot.outlier_mode.name}"] = plot.mplobj
        # combine figures
        utils.combine_figures(
            children,
            nsamples=CONDITION.nsample,
            outlier_modes=CONDITION.outlier_mode,
            output=self.filepath if save else None,
        )


@dataclass
class Paper:
    nrep: int = field(default=Default.NREP)
    seed: int = field(default=Default.SEED)
    figures: list[Figure] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        for simul_data in CONDITION.simul_data:
            for disp_type in CONDITION.disp_type:
                for frac_up in CONDITION.frac_up:
                    page = Figure(
                        simul_data=simul_data,
                        disp_type=disp_type,
                        frac_up=frac_up,
                        nrep=self.nrep,
                        seed=self.seed,
                    )
                    self.figures.append(page)

    def generate_datasets(self) -> list[Figure]:
        for figure in self.figures:
            figure.generate_datasets()
        return self.figures

    def calc_metrics(self) -> list[Figure]:
        for figure in self.figures:
            figure.calc_metrics()
        return self.figures

    @property
    def res2d(self) -> pd.DataFrame:
        return pd.concat([result.res2d for result in self.figures])

    def draw(self, save: bool = True) -> list[Figure]:
        for figure in self.figures:
            figure.draw(save=save)
        return self.figures

    def make(self, save: bool = True) -> list[Figure]:
        for figure in self.figures:
            figure.make(save=save)
        return self.figures
