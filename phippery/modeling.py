r"""
=================
Modeling
=================

A collection of interfaces for modeling binding enrichments,
and estimating significance.
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
from phippery.zscore import zscore_pids_binning
from phippery.zscore import compute_zscore

# TODO K: could you put the model equations in the docstrings?

def gamma_poisson_model(
    ds,
    starting_alpha=0.8,
    starting_beta=0.1,
    trim_percentile=99.9,
    data_table="size_factors",
    inplace=True,
    new_table_name="gamma_poisson_mlxp",
):
    """Fit a gamma poisson distribution and estimate the
    :math:`-log_{10}(p)` value, or *mlxp*
    for each sample-peptide enrichment in the dataset provided

    Note
    ----
    much of this source code is derived from 
    https://github.com/lasersonlab/phip-stat
    and written by Uri Laserson.


    Parameters
    ----------
    ds : xarray.DataSet
        The dataset you would like to fit to

    starting_alpha : float
        TODO K: decription

    starting_beta : float

    trim_percentile : float

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
        
    Returns
    -------
    tuple :
        the alpha, and beta values of the fit. If inplace is false, a copy of the new dataset
        is returned first.
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
    trim_percentile=100.,
    outlier_reject_scale=10.,
    data_table="size_factors",
    inplace=True,
    new_table_name="neg_binom_mlxp",
):
    """Fit a negative binomial distribution and estimate the
    :math:`-log_{10}(p)` value, or *mlxp*
    for each sample-peptide enrichment in the dataset provided.

    Parameters
    ----------
    ds : xarray.DataSet
        The dataset containing samples to estimate significance on.

    beads_ds : xarray.DataSet
        The dataset containing beads only control samples
        for which the distribution will be fit.

    nb_p : int
        TODO K: decription

    starting_beta : float

    outlier_reject_scale : float

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray

    Returns
    -------
    tuple :
        TODO K: description
    """
    #'nb_p' determines the relationship between mean and variance. Valid values
    #are 1 and 2 (sometimes called Type-1 and Type-2 Negative Binominal, respectively)

    #If 'inplace' parameter is True, then this function
    #appends a dataArray to ds which is indexed with the same coordinate dimensions as
    #'data_table'. If False, a copy of ds is returned with the appended dataArray

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
        (mu, alpha, var, size, prob) = fit_neg_binom(
            trimmed_data[i].compressed(), nb_p, outlier_reject_scale
        )
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


def zscore(
    ds,
    beads_ds,
    data_table='cpm',
    min_Npeptides_per_bin=300,  
    lower_quantile_limit=0.05,  
    upper_quantile_limit=0.95,  
    inplace=True,
    new_table_name='zscore'
):
    """Fit a negative binomial distribution and estimate the
    :math:`-log_{10}(p)` value, or *mlxp*
    for each sample-peptide enrichment in the dataset provided.

    Parameters
    ----------
    ds : xarray.DataSet
        The dataset containing samples to estimate significance on.

    beads_ds : xarray.DataSet
        for which the distribution will be fit.

    min_Npeptides_per_bin : int
        Mininum number of peptides per bin.

    lower_quantile_limit : float
        Counts below this quantile are ignored for computing mean, stddev.

    upper_quantile_limit : float
        Counts above this quantile are igonred for computing mean, stddev.

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray

    Returns
    -------
    tuple :
        TODO K: description
    """
    #This is a wrapper function for our xarray dataset.

    #for each sample in the dataset provided, compute Z-score following the method described
    #in the supplement to DOI:10.1126/science.aay6485

    #If 'inplace' parameter is True, then this function
    #appends a dataArray to ds which is indexed with the same coordinate dimensions as
    #'data_table'. If False, a copy of ds is returned with the appended dataArray
    #"""

    binning = zscore_pids_binning(beads_ds, data_table, min_Npeptides_per_bin)
    
    zscore_table = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    zs_df, mu_df, sigma_df = compute_zscore(ds, data_table, binning, lower_quantile_limit, upper_quantile_limit)
    zscore_table.loc[:, :] = zs_df

    if inplace:
        ds[new_table_name] = xr.DataArray(zscore_table)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(zscore_table)
        return ds_copy
