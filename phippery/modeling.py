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

    # TODO check that data of choice is in ds
    # TODO append to sample table the alpha and beta rates

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
