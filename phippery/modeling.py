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
from phippery.zscore import zscore_pids_binning
from phippery.zscore import compute_zscore


def gamma_poisson_model(
    ds,
    starting_alpha=0.8,
    starting_beta=0.1,
    trim_percentile=99.9,
    data_table="size_factors",
    inplace=True,
    new_table_name="gamma_poisson_mlxp",
):
    r"""Fit a Gamma distribution to determine Poisson rates
    per peptide for the non-specific binding background and estimate the
    :math:`-\log_{10}(p)` value, or *mlxp*,
    for each sample-peptide enrichment in the dataset provided.
    We use the following parameterization of the Gamma distribution:

    .. math::
        f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}

    The fit is performed on the distribution of peptide average counts to
    obtain :math:`\alpha` and :math:`\beta`. If there are :math:`n`
    samples involved in the fit, and for a given peptide with counts
    :math:`x_1, x_2, \ldots, x_n`, the background Poisson distribution
    is determined by the rate,

    .. math::
        \lambda = \frac{\alpha + \sum^n_{k=1} x_k}{\beta + n}

    Note
    ----
    much of this source code is derived from
    https://github.com/lasersonlab/phip-stat
    and written by Uri Laserson.


    Parameters
    ----------
    ds : xarray.DataSet
        The dataset you would like to fit to.

    starting_alpha : float
        Initial value for the shape parameter of the Gamma distribution.

    starting_beta : float
        Initial value for the rate parameter of the Gamma distribution.

    trim_percentile : float
        The percentile cutoff for removing peptides with very high counts.
        (e.g. a value of 98 means peptides in the highest 2% in counts 
        would be removed from the fit)
        This parameter is used to remove potential signal peptides that
        would bias the fit.

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


def zscore(
    ds,
    beads_ds,
    data_table="cpm",
    min_Npeptides_per_bin=300,
    lower_quantile_limit=0.05,
    upper_quantile_limit=0.95,
    inplace=True,
    new_table_name="zscore",
):
    r"""Calculate a Z-score of empirical enrichment relative to
    expected background mean CPM (:math:`\mu`) and stddev CPM (:math:`\sigma`)
    from beads-only samples,
    for each sample-peptide enrichment in the dataset provided.
    For a peptide with CPM :math:`n`, the Z-score is,

    .. math::
        z = \frac{n - \mu}{\sigma}

    Note
    ----
    This implementation follows the method described in the
    supplement to DOI:10.1126/science.aay6485.

    Parameters
    ----------
    ds : xarray.DataSet
        The dataset containing samples to estimate significance on.

    beads_ds : xarray.DataSet
        The dataset containing beads only control samples to estimate
        background means and stddevs.

    min_Npeptides_per_bin : int
        Mininum number of peptides per bin.

    lower_quantile_limit : float
        Counts below this quantile are ignored for computing background mean and stddev.

    upper_quantile_limit : float
        Counts above this quantile are igonred for computing background mean and stddev.

    data_table : str
        The name of the enrichment layer from which you would like to compute Z-scores.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray

    Returns
    -------
    None, xarray.DataSet :
        If inplace is false, a copy of the new dataset is returned.
    """
    # This is a wrapper function for our xarray dataset.

    # for each sample in the dataset provided, compute Z-score following the method described
    # in the supplement to DOI:10.1126/science.aay6485

    # If 'inplace' parameter is True, then this function
    # appends a dataArray to ds which is indexed with the same coordinate dimensions as
    #'data_table'. If False, a copy of ds is returned with the appended dataArray
    # """

    binning = zscore_pids_binning(beads_ds, data_table, min_Npeptides_per_bin)

    zscore_table = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    zs_df, mu_df, sigma_df = compute_zscore(
        ds, data_table, binning, lower_quantile_limit, upper_quantile_limit
    )
    zscore_table.loc[:, :] = zs_df

    if inplace:
        ds[new_table_name] = xr.DataArray(zscore_table)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(zscore_table)
        return ds_copy
