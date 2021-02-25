"""
@File: svd.py

@Author: Jared Galloway

This file contains a collection of functions which
uses Singular Value Decomposition to construct
Eigenassays and Eigenpeptides

Each samples can be projected into the eigenspace
and each peptide from the original data can then be projected
into the eigepeptide space.
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
from phippery.utils import id_coordinate_subset
from phippery.tidy import tidy_ds


def svd_aa_loc(
    ds,
    rank=1,
    data_table="enrichment",
    scaled_by_wt=False,
    protein_name_column="Protein",
    location_col="Loc",
    aa_sub_col="aa_sub",
    inplace=True,
    new_table_name="svd_rr",
):
    """
    compute singular value decomposition rank reduction
    on the aa / loc matrix by pivoting before computing decomposiion
    and re-shaping to add to the dataset.

    :param: r <int> Number of ranks in re-composition estimate.
    """

    low_rank_dt = copy.deepcopy(ds[data_table].to_pandas())

    for prot, prot_ds in iter_peptide_groups(ds, protein_name_column):

        for sid in prot_ds.sample_id.values:

            # grab the single sample ds
            rep_ds = prot_ds.loc[
                dict(
                    sample_id=[sid],
                    sample_metadata=["sample_ID"],
                    peptide_metadata=[aa_sub_col, location_col],
                )
            ]

            # melt
            tidy = tidy_ds(rep_ds)

            # Pivot so that we get the (aa X Loc)
            piv = tidy.pivot(index=aa_sub_col, columns=location_col, values=data_table)

            # Preserve the indices for population of new enrichment table
            piv_index = tidy.pivot(
                index=aa_sub_col, columns=location_col, values="peptide_id"
            )

            # compute rank reduction decompisition matrices
            U, S, V = svd(piv)

            # Grab the first X outer products in the finite summation of rank layers.
            low_rank = U[:, :rank] @ np.diag(S[:rank]) @ V[:rank, :]

            low_rank_piv = pd.DataFrame(low_rank, index=piv.index, columns=piv.columns)
            melted_values = pd.melt(low_rank_piv.reset_index(), id_vars=[aa_sub_col])
            melted_index = pd.melt(piv_index.reset_index(), id_vars=[aa_sub_col])
            melted_values["peptide_id"] = melted_index["value"]
            low_rank_dt.loc[melted_values["peptide_id"], sid] = melted_values[
                "value"
            ].values

    svd_rr_approx = xr.DataArray(low_rank_dt, dims=ds.counts.dims)

    if inplace:
        ds[new_table_name] = svd_rr_approx
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = svd_rr_approx
        return ds_copy


def compute_svd_decomposition(
    ds,
    data_table="counts",
    compute_assay_projections=False,
    compute_assay_correlations=False,
    # compute_peptide_projections=False,
    # compute_peptide_correlations=False,
):

    a = ds[f"{data_table}"].values
    U, S, V = svd(a, full_matrices=False)
    n_eigenvectors = len(ds.sample_id)

    z_jk = np.zeros([n_eigenvectors, n_eigenvectors])
    p_jk = np.zeros([n_eigenvectors, n_eigenvectors])
    for j in range(n_eigenvectors):
        for k in range(n_eigenvectors):
            z_jk[j, k] = np.dot(a[:, j], U[:, k])
            p_jk[j, k] = st.pearsonr(a[:, j], U[:, k])[0]

    sam_meta = ds.sample_table.to_pandas()
    for k in range(n_eigenvectors):
        sam_meta[f"Eigenassay-{k}-projection"] = z_jk[:, k]
        sam_meta[f"Eigenassay-{k}-correlation"] = p_jk[:, k]

    return (U, S, V), sam_meta
