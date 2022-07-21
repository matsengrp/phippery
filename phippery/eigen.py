r"""
=================
Eigen
=================

Defines a collection of functions which
uses Singular Value Decomposition to construct
Eigenassays and Eigenpeptides

Specifically, 
"""


import numpy as np
from numpy.linalg import svd
import xarray as xr
import pandas as pd
import scipy.stats as st
import itertools
import copy

from phippery.utils import iter_peptide_groups
from phippery.utils import iter_sample_groups


def eigenassay_projections(
    ds,
    data_table="counts",
    compute_correlations=False,
    return_raw_decomposition=False,
    return_eigenassay_meta=False,
    n_eigenvectors=None,
):

    """Compute the Singular Value Decomposition
    of the enrichment data before projecting each
    sample into the first n eigenvector ("eigenassay") 
    dimensions in the dataset.

    Parameters
    ----------
    ds : xarray.DataSet
        The dataset you would like to perform eigen decomposition on

    Returns
    -------
    dict :
        1. the raw "economy" SVD decomposition matrices
        2. the eigenassays tied with the peptide metadata included in `ds`

    Example
    -------
        TODO
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    n_eigenvectors = len(ds.sample_id) if n_eigenvectors is None else n_eigenvectors
    if n_eigenvectors > len(ds.sample_id):
        raise ValueError(
            f"n_eigenvectors must be less than or equal to the number of samples in the dataset"
        )

    a = ds[f"{data_table}"].values
    U, S, V = svd(a, full_matrices=False)
    sam_meta = ds.sample_table.to_pandas()

    z_jk = np.zeros([n_eigenvectors, n_eigenvectors])
    if compute_correlations:
        p_jk = np.zeros([n_eigenvectors, n_eigenvectors])
    for j in range(n_eigenvectors):
        for k in range(n_eigenvectors):
            z_jk[j, k] = np.dot(a[:, j], U[:, k])
            if compute_correlations:
                p_jk[j, k] = st.pearsonr(a[:, j], U[:, k])[0]

    for k in range(n_eigenvectors):
        sam_meta[f"Eigenassay-{k}-projection"] = z_jk[:, k]
        if compute_correlations:
            sam_meta[f"Eigenassay-{k}-correlation"] = p_jk[:, k]

    ret = {"sample_eigenassay_projections": sam_meta}
    if return_raw_decomposition:
        ret["raw_decomposition"] = (U, S, V)
    if return_eigenassay_meta:
        p_table = ds.peptide_table.to_pandas()
        for rank in range(n_eigenvectors):
            p_table[f"Column_Eigenassay_{rank}"] = U[:, rank].flatten()
        ret["eigenassay_meta"] = p_table

    return ret
