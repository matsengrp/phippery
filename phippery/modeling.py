"""
@File: modeling.py

@Author: Jared Galloway
"""

import numpy as np
import xarray as xr
import pandas as pd
import itertools
import copy
import scipy.stats as st
from phippery.gampois import fit_gamma
from phippery.gampois import gamma_poisson_posterior_rates
from phippery.gampois import mlxp_gamma_poisson
from phippery.negbinom import fit_neg_binom
from phippery.negbinom import mlxp_neg_binom


def gamma_poisson_model(
    ds,
    starting_alpha=0.8,
    starting_beta=0.1,
    trim_percentile=99.9,
    data_table="size_factors",
    inplace=True,
    new_table_name="gamma_poisson_mlxp",
):
    """
    This is a wrapper function for our xarray dataset.
    The original code can be found here:
    https://github.com/lasersonlab/phip-stat
    written by Uri Laserson.

    for each sample in the dataset provided,
    compute -log10(pval) for counts matrix
    counts is DataFrame; assumed columns are normalized to some size factor.

    The function will return a tuple containing the alpha, and beta values of the
    fit.

    If 'inplace' parameter is True, then this function
    appends a dataArray to ds which is indexed with the same coordinate dimensions as
    'data_table'. If False, a copy of ds is returned with the appended dataArray
    """

    print(
        f"Running new modeling with starting alpha = {starting_alpha}, beta = {starting_beta}"
    )

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    counts = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    counts = counts.round(2)

    upper_bound = st.scoreatpercentile(counts.values, trim_percentile)
    trimmed_means = np.ma.mean(
        np.ma.masked_greater(counts.values, upper_bound), axis=1
    ).data
    alpha, beta = fit_gamma(
        trimmed_means, starting_alpha=starting_alpha, starting_beta=starting_beta
    )
    background_rates = gamma_poisson_posterior_rates(counts, alpha, beta, upper_bound)
    counts.loc[:, :] = mlxp_gamma_poisson(counts, background_rates)

    if inplace:
        ds[new_table_name] = xr.DataArray(counts)
        return (alpha, beta)
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(counts)
        return (alpha, beta), ds_copy


def neg_binom_model(
    ds,
    beads_ds,
    nb_p=2,
    trim_percentile=99.9,
    data_table="size_factors",
    inplace=True,
    new_table_name="neg_binom_mlxp",
):
    """
    This is a wrapper function for our xarray dataset.

    for each sample in the dataset provided,
    compute -log10(pval) for counts matrix
    counts is DataFrame; assumed columns are normalized to some size factor.

    'nb_p' determines the relationship between mean and variance. Valid values
    are 1 and 2 (sometimes called Type-1 and Type-2 Negative Binominal, respectively)

    If 'inplace' parameter is True, then this function
    appends a dataArray to ds which is indexed with the same coordinate dimensions as
    'data_table'. If False, a copy of ds is returned with the appended dataArray
    """

    # TODO check that data of choice is in ds

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    beads_counts = copy.deepcopy(beads_ds[f"{data_table}"].to_pandas())
    upper_bound = st.scoreatpercentile(beads_counts.values, trim_percentile)
    trimmed_data = np.ma.masked_greater(beads_counts.values, upper_bound)

    nb_mu = []
    nb_alpha = []
    nb_var = []
    nb_size = []
    nb_prob = []
    for i in range(beads_counts.shape[0]):
        (mu, alpha, var, size, prob) = fit_neg_binom(trimmed_data[i].compressed(), nb_p)
        nb_mu.append(mu)
        nb_alpha.append(alpha)
        nb_var.append(var)
        nb_size.append(size)
        nb_prob.append(prob)

    counts = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    counts = counts.round(2)
    counts.loc[:, :] = mlxp_neg_binom(counts, nb_size, nb_prob)

    if inplace:
        ds[new_table_name] = xr.DataArray(counts)
        return (nb_size, nb_prob)
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(counts)
        return (nb_size, nb_prob), ds_copy
