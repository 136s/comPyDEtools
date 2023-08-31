#!/usr/bin/env python

"""
useful functions and classes
"""

from copy import deepcopy
import platform
import subprocess
from typing import Generator

from matplotlib.axes import Axes
from matplotlib.figure import Figure as mplFigure
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.metrics import (
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
    f1_score,
    cohen_kappa_score,
)
import skunk

from .const import Default, Metrics, MetricsInput, Outlier, StrPath
from .generation import GENE_DESCRIPTION_COLNAME, UP_DEG, DN_DEG


FONT: dict[str, dict[str, list[str]]] = {
    "Windows": ["Arial", "Nimbus Sans", "Helvetica"],
    "Linux": ["Nimbus Sans", "Arial", "Helvetica"],
    "Darwin": ["Helvetica", "Nimbus Sans", "Arial"],
}


def run_commands(cmds: list) -> Generator[str, None, None]:
    for cmd in cmds:
        p = subprocess.run(cmd, shell=True)
        yield p.stdout


def load_metrics_input(
    de_output_path: StrPath,
    de_true_colname: str = GENE_DESCRIPTION_COLNAME,
    de_score_colname: str = "fdr",
    de_score_threshold: float = 0.1,
) -> MetricsInput:
    metrics_input = (
        pd.read_csv(de_output_path)
        .rename(
            columns={
                de_true_colname: "y_true",
                de_score_colname: "y_score",
            }
        )[["y_true", "y_score"]]
        .dropna()
        .to_dict(orient="list")
    )
    metrics_input["y_true"] = [r in (UP_DEG, DN_DEG) for r in metrics_input["y_true"]]
    metrics_input["y_score"] = [1 - s for s in metrics_input["y_score"]]
    metrics_input["y_pred"] = [
        s > 1 - de_score_threshold for s in metrics_input["y_score"]
    ]
    return metrics_input


def calc_metrics(
    metrics_type: Metrics,
    y_true: list[float] | None = None,
    y_score: list[float] | None = None,
    y_pred: list[float] | None = None,
) -> float:
    match metrics_type:
        case Metrics.auc:
            metrics_value = roc_auc_score(y_true=y_true, y_score=y_score)
        case Metrics.tpr:  # = TP/(TP+FN)
            metrics_value = recall_score(
                y_true=y_true, y_pred=y_pred, labels=[True, False]
            )
        case Metrics.fdr:  # = FP/(TP+FP)
            metrics_value = 1 - precision_score(
                y_true=y_true, y_pred=y_pred, zero_division=0, labels=[True, False]
            )
        case Metrics.cutoff:
            _fpr, _tpr, _thresholds = roc_curve(
                y_true=y_true, y_score=y_score, pos_label=True
            )
            metrics_value = 1 - _thresholds[(_fpr**2 + (1 - _tpr) ** 2).argmin()]
        case Metrics.f1score:
            metrics_value = f1_score(y_true=y_true, y_pred=y_pred, labels=[True, False])
        case Metrics.kappa:
            metrics_value = cohen_kappa_score(
                y1=y_true, y2=y_pred, labels=[True, False]
            )
        case _:
            raise ValueError
    return metrics_value


def set_font() -> None:
    plt.rcParams["svg.fonttype"] = "none"
    default_family = "font." + plt.rcParams["font.family"][0]
    plt.rcParams[default_family] = (
        FONT[platform.system()] + plt.rcParams[default_family]
    )


def draw_plot(
    data: pd.DataFrame,
    plot_title: str | None = None,
    height: float = 2,
    palette: list[str] = Default.PALETTE,
    output: StrPath | None = None,
) -> sns.FacetGrid:
    aspect = len(data["method"].unique()) / 5
    sns.set_style("darkgrid", rc={"axes.facecolor": "ebebeb"})
    set_font()
    g: sns.FacetGrid = sns.catplot(
        data,
        x="method",
        y="value",
        row="metrics",
        col="pde",
        height=height,
        aspect=aspect,
        kind="box",
        palette=palette,
        sharex="col",
        sharey="row",
        showcaps=False,
        flierprops={"marker": "x"},
        linewidth=height / 2,
        margin_titles=True,
    )
    if plot_title:
        g.figure.suptitle(plot_title)
    g.set_axis_labels("", "")
    g.set_titles(col_template="pDE = {col_name}%", row_template="{row_name}", size=14)
    for ax in g.axes.flat:
        ax.grid(which="both", axis="both", color="white", linewidth=2)
        plt.setp(ax.get_xticklabels(), rotation=90)
    g.tight_layout()
    if output:
        g.savefig(output)
    plt.close()
    return g


def combine_figures(
    plots: dict[str, mplFigure],
    nsamples: list[int],
    outlier_modes: list[Outlier],
    output: StrPath = None,
) -> str:
    # get figure size
    w, h = 0, 0
    for plot in plots.values():
        w = max(w, plot.get_tightbbox().width)
        h = max(h, plot.get_tightbbox().height)
    # prepare parent figure
    ncols = len(nsamples)
    nrows = len(outlier_modes)
    fig, axs = plt.subplots(nrows, ncols, figsize=(w * ncols, h * nrows))
    fig.tight_layout(pad=0)
    # set skunk gid
    for x, col in enumerate(nsamples):
        for y, row in enumerate(outlier_modes):
            if ncols == 1 and nrows == 1:
                ax: Axes = axs
            elif ncols == 1:
                ax: Axes = axs[y]
            elif nrows == 1:
                ax: Axes = axs[x]
            else:
                ax: Axes = axs[y][x]
            ax.axis("off")
            skunk.connect(ax, f"{col}spc_{row.name}")
    # insert mpl figures to parent figure
    combined = skunk.insert(deepcopy(plots))
    if output:
        with open(output, "w") as f:
            f.write(combined)
    plt.close()
    return combined
