r"""
=================
Eigen
=================

A method for MDS via Singular Value Decomposition 
on phip-seq datasets.
"""

import numpy as np
from numpy.linalg import svd
import xarray as xr
import pandas as pd
import scipy.stats as st
import itertools
import copy

def eigenassay_projections(
    ds,
    data_table="counts",
    compute_correlations=False,
    return_raw_decomposition=False,
    n_eigenvectors=None,
):

    r"""Compute the Singular Value Decomposition
    of the enrichment data before projecting each
    sample into the first n eigenvectors ("eigenassays")
    in the dataset.

    Concretely, given a Matrix of, :math:`X` enrichments in the
    `phippery` dataset with shape (peptides, samples). We compute
    the decomposition :math:`X = USV^{T}`

    The principal axes in feature space are then represented by the columns of 
    :math:`V` and represent the direction of maximum variance in the data. 
    The sample projections into this space are then computed and tied to the
    sample annotation in the returned dictionary.

    Parameters
    ----------
    ds : xarray.DataSet
        The dataset you would like to perform eigen-decomposition on.

    data_table : str
        The name of the enrichment layer in the phippery dataset to perform
        the operation on.

    compute_correlations : bool
        If true, compute the correlation of a sample's true binding enrichments to
        the n'th principal axes. These will be added in the same fashion as the sample
        scores that are appended to the sample table.

    return_raw_decomposition : bool
        If true, include the raw :math:`USV^{T}` decomposition of the enrichment
        matrix specified

    n_eigenvectors : int
        the number of projections "eigenassay dimensions" to include.


    Returns
    -------
    dict :
        1. The eigenassay projects tied with the appended to sample annotations
            included in `ds`.
        2. (optional) The raw "economy" SVD decomposition matrices.

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

    return ret
