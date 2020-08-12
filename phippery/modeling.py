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


def gamma_poisson_model(ds, trim_percentile=99.9, data_table="size_factors"):
    """
    This is a wrapper function for out xarray dataset.
    The original code can be found here:
    https://github.com/lasersonlab/phip-stat
    written by Uri Laserson.

    for each sample in the dataset provided,
    compute -log10(pval) for counts matrix
    counts is DataFrame; assumed columns are normalized to some size factor.

    The function will return another xarray dataset where another DataArray
    has been added to the dataset named `mlxp` that is indexed with the
    same dimentional coordinates (peptide_id, sample_id) as the counts table.
    """

    # TODO check that data of choice is in ds
    # TODO append to sample table the alpha and beta rates

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    counts = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    counts = counts.round(2)
    # ret_sample_table = ds.sample_table.to_pandas()
    # for new_column in ["gp_alpha", "gp_beta", "gp_upper_bound"]:
    #    assert new_column not in ret_sample_table.columns
    #    ret_sample_table[new_column] = np.zeros(len(ret_sample_table))

    upper_bound = st.scoreatpercentile(counts.values, trim_percentile)
    trimmed_means = np.ma.mean(
        np.ma.masked_greater(counts.values, upper_bound), axis=1
    ).data
    alpha, beta = fit_gamma(trimmed_means)
    # ret_sample_table.loc[sample_id, "gp_alpha"] = alpha
    # ret_sample_table.loc[sample_id, "gp_beta"] = beta
    background_rates = gamma_poisson_posterior_rates(counts, alpha, beta, upper_bound)
    counts.loc[:, :] = mlxp_gamma_poisson(counts, background_rates)
    ds[f"gamma_poisson_mlxp_{data_table}"] = xr.DataArray(counts)
