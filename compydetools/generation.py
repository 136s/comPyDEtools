#!/usr/bin/env python

"""
Porting R's `compareDEtools::SyntheticDataSimulation()` function.
"""

from copy import deepcopy
from importlib import resources
from typing import Literal

import numpy as np
import pandas as pd

from . import data
from .const import Default, Simul, Disp, Outlier


# dictionaly of simulation parameters
DATASET_PARAMETERS: dict[
    Literal[
        "k_total_mean",
        "k_total_disp",
        "k_cancer_mean",
        "k_cancer_disp",
        "k_normal_mean",
        "k_normal_disp",
        "b_total_mean",
        "b_total_disp",
        "b_C_mean",
        "b_C_disp",
        "b_D_mean",
        "b_D_disp",
    ],
    pd.Series,
]


def order(data: list[int] | pd.Series) -> list[int]:
    """returns a permutation which rearranges data into ascending or descending order

    Args:
        data (list[int] | pd.Series): Data list.

    Returns:
        list[int]: List of indices in ascending order.

    Notes:
        Porting R's `base::order()` function.
    """
    return sorted(range(len(data)), key=lambda i: data[i])


def rnbinom(
    mean: float,
    disp: float,
    size: int | tuple[int] | None = None,
    rng: np.random.Generator | None = None,
) -> float | np.ndarray:
    """Get random generation for the negative binomial distribution with parameters

    Args:
        mean (float): mean. `mu` in R.
        disp (float): dispersion. `1/size` in R.
        size (int | tuple[int] | None, optional): number of observations. `n` in R.
            Defaults to None.
        rng (np.random.Generator | None, optional): random.Generator object.
            Defaults to None.

    Returns:
        float | np.ndarray: Samples drawn from a parameterized negative binomial
            distribution, where each sample represents the number of failures that
            occurred before achieving `1/disp` successes.

    Notes:
        Created with compatibility to R's `stats::rnbinom(n, size, prob, mu)`.
    """
    if rng is None:
        rng = np.random.default_rng()
    return rng.negative_binomial(n=1 / disp, p=1 / (mean * disp + 1), size=size)


def get_disp(
    mean: float,
    mean_condition: pd.Series,
    disp_condition: pd.Series,
    rng: np.random.Generator | None = None,
) -> float:
    """
    Get random dispersion value based on the read count size (mean).

    Args:
        mean (float): Read count size (mean).
        mean_condition (pd.Series): Read count vector of a certain sample condition.
        disp_condition (pd.Series): Dispersion vector of a certain sample condition.
        rng (np.random.Generator | None, optional): random.Generator object.
            Defaults to None.

    Returns:
        float: Random dispersion value.

    Notes:
        Porting R's `compareDEtools::getDisp()` function.
    """
    # Get disp for elements where mean is within ±20.
    pool = disp_condition[
        (mean_condition > (mean - 20)) & (mean_condition < (mean + 20))
    ]
    if len(pool) == 0:
        # If no match is found, retrieve disp for the element with the closest mean
        # (choose only one if there are multiple means with the same closeness).
        value = disp_condition[np.argmin(np.abs(mean_condition - mean))]
    else:
        if rng is None:
            rng = np.random.default_rng()
        # If there are one or more matches, choose randomly.
        value = rng.choice(pool, 1)[0]
    return value


def generate_dataset_parameter() -> None:
    """Get a list containing estimated mean and dispersion parameters

    Notes:
        The parameters are precomputed values obtained using
        `data/generateDatasetParameter/SimulationDataGeneration.R`.
    """
    global DATASET_PARAMETERS
    try:
        DATASET_PARAMETERS
    except NameError:
        k_params_path = resources.files(data).joinpath("k_params.csv")
        b_params_path = resources.files(data).joinpath("b_params.csv")
        DATASET_PARAMETERS = dict(
            **pd.read_csv(k_params_path, index_col=0).to_dict(orient="series"),
            **pd.read_csv(b_params_path, index_col=0).to_dict(orient="series"),
        )


def synthetic_data_simulation(
    simul_data: Simul = Simul.first,
    random_sampling: bool = True,
    large_sample: bool = False,
    fixedfold: bool = False,
    nsample: int = 3,
    ngenes: int = 5000,
    nde: int = 250,
    frac_up: float = 0.5,
    disp_type: Disp = Disp.first,
    outlier_mode: Outlier = Outlier.first,
    ro_prop: int = 5,
    seed: int = Default.SEED,
) -> pd.DataFrame:
    """
    Generate synthetic count data for analysis.

    Args:
        simul_data (Simul, optional): dataset will be used for simulation data
            generation. Simul("KIRC") for KIRC dataset. Simul("Bottomly") for Bottomly
            dataset. Simul("mKdB") for hybrid dataset combining mean of KIRC and
            dispersion of Bottomly datatset. Simul("mBdK") for hybrid dataset combining
            mean of Bottomly and dispersion of KIRC datatset. Defaults to Simul.first.
        random_sampling (bool, optional): If True, introduce randomness among samples
            within the group. Defaults to True.
        large_sample (bool, optional): If True, increase the count by 5-10 times for
            only 3% of the non-outlier samples (OL3), only when outlier_mode == 'OS'.
            Defaults to False.
        fixedfold (bool, optional): whether this dataset is generated by fold changes
            with fixed values or random folds following exponential distribution, only
            when simul_data == 'KIRC'. Defaults to False.
        nsample (int, optional): number of samples for each sample group. Defaults to 3.
        ngenes (int, optional): the number of total gene in the synthetic data. Defaults
            to 5000.
        nde (int, optional): number of generated DE genes in the synthetic data.
            Defaults to 250.
        frac_up (float, optional): Proportion of upregulated DE genes in the synthetic
            data.  Defaults to 0.5.
        disp_type (Disp, optional): how is the dispersion parameter assumed to be for
            each condition to make a synthetic data. Defaults to Disp.first.
        outlier_mode (Outlier, optional): outlier mode for simulation data generation.
            "D" for basic simulation (not adding outliers). "R" for adding 5% of random
            outlier. "OS" for adding outlier sample to each sample group. "DL" for
            decreasing KIRC simulation dispersion 22.5 times (similar to SEQC data
            dispersion) to compare with SEQC data. Defaults to Outlier.first.
        ro_prop (int, optional): random outlier proportion percentage that we generate
            dataset with, only when outlier_mode == 'R'. Defaults to 5.
        seed (int, optional): seed of random number generation.
            Defaults to Default.SEED.

    Returns:
        pd.DataFrame: synthetic count
    """
    # relation
    # control : k_normal : b_D : suffix = 1
    # treatment : k_cancer : b_C : suffix = 2

    # set Random Generator
    rng = np.random.default_rng(seed)

    # Check the mode of outliers
    if not isinstance(outlier_mode, Outlier):
        raise ValueError(
            'outlier_mode must be "D" (DE variation test), '
            '"R" (random outlier test), '
            '"OS" (dispersion outlier sample test), or '
            '"DL" (dispersion lowered test)'
        )

    # set dataset parameters
    generate_dataset_parameter()
    global DATASET_PARAMETERS
    k_normal_disp: pd.Series = DATASET_PARAMETERS["k_normal_disp"].copy(deep=True)
    k_cancer_disp: pd.Series = DATASET_PARAMETERS["k_cancer_disp"].copy(deep=True)
    b_D_disp: pd.Series = DATASET_PARAMETERS["b_D_disp"].copy(deep=True)
    b_C_disp: pd.Series = DATASET_PARAMETERS["b_C_disp"].copy(deep=True)
    k_normal_mean: pd.Series = DATASET_PARAMETERS["k_normal_mean"].copy(deep=True)
    k_cancer_mean: pd.Series = DATASET_PARAMETERS["k_cancer_mean"].copy(deep=True)
    b_D_mean: pd.Series = DATASET_PARAMETERS["b_D_mean"].copy(deep=True)
    b_C_mean: pd.Series = DATASET_PARAMETERS["b_C_mean"].copy(deep=True)
    k_total_mean: pd.Series = DATASET_PARAMETERS["k_total_mean"].copy(deep=True)
    k_total_disp: pd.Series = DATASET_PARAMETERS["k_total_disp"].copy(deep=True)
    b_total_mean: pd.Series = DATASET_PARAMETERS["b_total_mean"].copy(deep=True)
    b_total_disp: pd.Series = DATASET_PARAMETERS["b_total_disp"].copy(deep=True)

    # by the simulation parameters, get the average counts of samples for each gene.
    if simul_data in (Simul.KIRC, Simul.mKdB):
        sample_mean1 = k_total_mean.sample(n=ngenes, random_state=rng)
        if simul_data == Simul.mKdB:
            sub_sample_mean1 = b_total_mean.sample(n=ngenes, random_state=rng)
            sub_sample_mean2 = deepcopy(sub_sample_mean1)
    elif simul_data in (Simul.Bottomly, Simul.mBdK):
        sample_mean1 = b_total_mean.sample(n=ngenes, random_state=rng)
        if simul_data == Simul.mBdK:
            sub_sample_mean1 = k_total_mean.sample(n=ngenes, random_state=rng)
            sub_sample_mean2 = deepcopy(sub_sample_mean1)
    sample_mean2 = deepcopy(sample_mean1)

    # set deg mean
    if nde != 0:
        if fixedfold:
            # set the ratio of up DEG to 2/3 to enhance comparability with SEQC
            if simul_data != Simul.KIRC:
                raise ValueError(
                    "Simulation with fixed fold must be based on KIRC dataset."
                )
            factor1 = 1.15  # up
            factor2 = 1.3  # up
            factor3 = 1.6  # dn
            one_third = round(nde / 3)
            num_updeg = round(2 * nde / 3)
            sample_mean2[:one_third] *= factor1
            sample_mean2[one_third:num_updeg] *= factor2
            sample_mean2[num_updeg:nde] /= factor3
        else:
            # generate up-regulated DEG ratio with the specified value
            # calculate indices of DEG
            num_updeg = round(nde * frac_up)
            upindex = np.arange(num_updeg)
            dnindex = np.arange(num_updeg, nde)
            # vary the fold change based on the sample count
            # and add random numbers from exponential distribution exp(1)
            factor1 = rng.exponential(size=len(upindex))
            factor2 = rng.exponential(size=len(dnindex))
            if nsample <= 3:  # 3 sample fold change 1.5~
                factor1 += 1.5
                factor2 += 1.5
            elif nsample <= 5:  # 5 sample fold change 1.3~
                factor1 += 1.3
                factor2 += 1.3
            else:  # 10 sample fold change 1.2~
                factor1 += 1.2
                factor2 += 1.2
            # calculate mean of DEG
            sample_mean2[upindex] *= factor1
            sample_mean2[dnindex] /= factor2
            if simul_data in (Simul.mBdK, Simul.mKdB):
                sub_sample_mean2[upindex] *= factor1
                sub_sample_mean2[dnindex] /= factor2

    # set the reference mean and dispersion for referencing
    if simul_data in (Simul.KIRC, Simul.mBdK):
        mean_condition1 = k_normal_mean
        mean_condition2 = k_cancer_mean
        disp_condition1 = k_normal_disp
        disp_condition2 = k_cancer_disp
        disp_total = k_total_disp
    elif simul_data in (Simul.Bottomly, Simul.mKdB):
        mean_condition1 = b_D_mean
        mean_condition2 = b_C_mean
        disp_condition1 = b_D_disp
        disp_condition2 = b_C_disp
        disp_total = b_total_disp

    # generate dispersion
    if disp_type == Disp.different:
        if simul_data in (Simul.KIRC, Simul.Bottomly):
            # when setting each disp for KIRC or Bottomly, randomly get disp based on
            # similar average expression from normal/cancer or D/C groups, respectively
            sample_disp1 = sample_mean1.apply(
                get_disp,
                mean_condition=mean_condition1,
                disp_condition=disp_condition1,
                rng=rng,
            )
            sample_disp2 = sample_mean1.apply(
                get_disp,
                mean_condition=mean_condition2,
                disp_condition=disp_condition2,
                rng=rng,
            )
        elif simul_data in (Simul.mKdB, Simul.mBdK):
            # setting each disp for mKdB or mBdK
            sample_disp1 = sub_sample_mean1.sort_values().apply(
                get_disp,
                mean_condition=mean_condition1,
                disp_condition=disp_condition1,
                rng=rng,
            )[order(order(sample_mean1))]
            sample_disp2 = sub_sample_mean2.sort_values().apply(
                get_disp,
                mean_condition=mean_condition2,
                disp_condition=disp_condition2,
                rng=rng,
            )[order(order(sample_mean2))]
    elif disp_type == Disp.same:
        # set same disp
        if simul_data in (Simul.KIRC, Simul.Bottomly):
            sample_disp1 = disp_total[sample_mean1.index]
        elif simul_data in (Simul.mKdB, Simul.mBdK):
            sample_disp1 = disp_total[sub_sample_mean1.index[order(sub_sample_mean1)]][
                order(order(sample_mean1))
            ]
        sample_disp2 = deepcopy(sample_disp1)

    # generate counts
    counts = np.zeros((ngenes, 2 * nsample))
    # calculate sample size
    one_third = round(nsample / 3)
    two_thirds = nsample - one_third
    four_thirds = nsample + one_third
    # Branching based on outlier sample mode
    if outlier_mode == Outlier.OS or large_sample:
        # increase disp by 5 times for 1/3 samples in each group (OL2 on paper)
        for i in range(ngenes):
            # outlier sample of treatment group
            counts[i, :one_third] = rnbinom(
                mean=sample_mean2[i], disp=5 * sample_disp2[i], size=one_third, rng=rng
            )
            # good sample of treatment group
            counts[i, one_third:nsample] = rnbinom(
                mean=sample_mean2[i], disp=sample_disp2[i], size=two_thirds, rng=rng
            )
            # outlier sample of control group
            counts[i, nsample:four_thirds] = rnbinom(
                mean=sample_mean1[i], disp=5 * sample_disp1[i], size=one_third, rng=rng
            )
            # good sample of control group
            counts[i, four_thirds:] = rnbinom(
                mean=sample_mean1[i], disp=sample_disp1[i], size=two_thirds, rng=rng
            )
        if large_sample:
            # increase counts by 5-10 times for only 3% of normal samples (OL3 on paper)
            # generate a set of uniform random numbers between 0 and 100
            ro = rng.uniform(low=0, high=100, size=ngenes * 2 * nsample)
            # find indices less than 3 to generate 3% outlier
            index_outlier = np.where(ro < 3)[0]
            # exact good sample indices
            index_outlier = index_outlier[
                np.where(
                    (  # treatment group
                        (index_outlier >= ngenes * one_third)
                        & (index_outlier < ngenes * nsample)
                    )
                    | (  # control group
                        (index_outlier >= ngenes * four_thirds)
                        & (index_outlier < ngenes * 2 * nsample)
                    )
                )
            ]
            # multiply counts for outlier by a uniform random number between 5 and 10
            counts_1d = counts.T.flatten()
            counts_1d[index_outlier] *= rng.uniform(
                low=5, high=10, size=len(index_outlier)
            )
            counts = np.round(counts_1d.reshape(counts.T.shape).T)
    elif outlier_mode == Outlier.DL:
        # to compare KIRC with SEQC, reduce disp by a factor of 1/22.5
        for i in range(ngenes):
            # treatment group
            counts[i, :nsample] = rnbinom(
                mean=sample_mean2[i],
                disp=sample_disp2[i] / 22.5,
                size=nsample,
                rng=rng,
            )
            # control group
            counts[i, nsample:] = rnbinom(
                mean=sample_mean1[i],
                disp=sample_disp1[i] / 22.5,
                size=nsample,
                rng=rng,
            )
    else:
        # not an outlier sample (OS; OL2), not a dispersion reduced sample (DL),
        # nor a large sample (OL3)
        if random_sampling:
            # if introducing randomness among samples within the group,
            # prepare a scaling factor
            rand1 = rng.uniform(low=0.7, high=1.3, size=nsample)
            rand2 = rng.uniform(low=0.7, high=1.3, size=nsample)
        else:
            rand1 = np.ones(nsample)
            rand2 = np.ones(nsample)
        for i in range(ngenes):
            # treatment group
            counts[i, :nsample] = [
                rnbinom(mean=sample_mean2[i] * scale, disp=sample_disp2[i], rng=rng)
                for scale in rand1
            ]
            # control group
            counts[i, nsample:] = [
                rnbinom(mean=sample_mean1[i] * scale, disp=sample_disp1[i], rng=rng)
                for scale in rand2
            ]

    # generate outlier
    if outlier_mode == Outlier.R:
        # increase counts of specified % by 5-10 times (OL1)
        ro = rng.uniform(low=0, high=100, size=ngenes * 2 * nsample)
        index_outlier = np.where(ro < ro_prop)[0]
        # multiply counts for outlier by a uniform random number between 5 and 10
        counts_1d = counts.T.flatten()
        counts_1d[index_outlier] *= rng.uniform(low=5, high=10, size=len(index_outlier))
        counts = np.round(counts_1d.reshape(counts.T.shape).T)

    # modify Gene IDs, sample names, Gene Symbols and Description
    count_matrix = pd.DataFrame(counts, dtype=int)
    count_matrix.index = count_matrix.index.to_series().add(1).rename("Gene_ID")
    count_matrix.columns = [f"TRT-{s+1}" for s in range(nsample)] + [
        f"CTRL-{s+1}" for s in range(nsample)
    ]
    gene_info = count_matrix.apply(lambda v: f"LOC{v.name}", axis=1).to_frame(
        "Gene_Symbol"
    )
    gene_info["Description"] = (
        ["up"] * num_updeg
        + ["dn"] * (nde - num_updeg)
        + ["ns"] * (len(gene_info) - nde)
    )

    return gene_info.join(count_matrix)
