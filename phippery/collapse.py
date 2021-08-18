"""
@File: collapse.py

@Author: Jared Galloway

This file contains code for collapsing a phip dataset
by sample groups
"""


# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as st

# built-in python3
from collections import defaultdict
import itertools

# local
from phippery.utils import iter_sample_groups


def throw_mismatched_features(df: pd.DataFrame, by: list):
    """
    When you collapse by some set of columns in the dataframe,
    keep only features which homogeneous within groups.

    This is similar to 'DataFrameGroupby.first()' method,
    but instead of keeping the first factor level appearing for each group
    category, we only throw any features which are note homogeneous within groups.
    You could achieve the same functionality by listing features you know to be
    homogeneous in the 'by' parameter.
    """

    # Find out which collapse features are shared within groups
    collapsed_sample_metadata = defaultdict(list)
    for i, (group, group_df) in enumerate(df.groupby(by)):
        for column, value in group_df.iteritems():
            v = value.values
            if np.all(v == v[0]) or np.all([n != n for n in v]):
                collapsed_sample_metadata[column].append(v[0])

    # Throw out features that are not shared between groups
    to_throw = [
        key for key, value in collapsed_sample_metadata.items() if len(value) < i + 1
    ]
    [collapsed_sample_metadata.pop(key) for key in to_throw]
    return pd.DataFrame(collapsed_sample_metadata)


def mean_pw_correlation_by_group(ds, by, data_tables="all"):
    """A wrapper for computing pw cc within groups defined
    with the 'by' parameter. Merges the correlations into a
    single table"""
    # NOTE TODO

    """
    # Compute mean pw cc on all possible data tables
    if data_tables == 'all':
        data_tables = list(
                set(ds.data_vars) - set(['sample_table', 'peptide_table'])
                )

    corr_tables = [
        pairwise_correlation_by_sample_group(ds, by, data_table)
        for data_table in data_tables
    ]

        #assert len(pw_cc) == len(csm)
        collapsed_anno_table = collapsed_anno_table.merge(
                pw_cc, left_index=True, right_index=True
                )
    """
    pass


def pairwise_correlation_by_group(ds, group, data_table="counts", column_prefix=None):

    """Computes pairwise cc for all
    sample in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of samples
    in the group."""

    # TODO, how could this

    if group not in ds.sample_metadata.values:
        raise ValueError(f"{group} does not exist in sample metadata")

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    data = ds[f"{data_table}"].to_pandas()

    groups, pw_cc, n = [], [], []
    for s_group, group_ds in iter_sample_groups(ds, group):
        groups.append(int(s_group))
        n.append(len(group_ds["sample_id"].values))

        if len(group_ds["sample_id"].values) < 2:
            pw_cc.append(1.0)
            continue
        correlations = []
        for sample_ids in itertools.combinations(group_ds["sample_id"].values, 2):
            sample_0_enrichment = data.loc[:, sample_ids[0]]
            sample_1_enrichment = data.loc[:, sample_ids[1]]
            correlation = (
                st.pearsonr(sample_0_enrichment, sample_1_enrichment)[0]
                if np.any(sample_0_enrichment != sample_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(round(sum(correlations) / len(correlations), 5))

    if column_prefix is None:
        column_prefix = f"{group}_{data_table}"

    ret = pd.DataFrame(
        {
            f"{group}": groups,
            f"{column_prefix}_pw_cc": np.array(pw_cc).astype(np.float64),
            f"{column_prefix}_n_reps": np.array(n).astype(np.int),
        }
    ).set_index(group)

    return ret


def collapse_groups(
    ds, by, collapse_dim="sample", agg_func=np.mean, compute_pw_cc=False
):
    """
    Collapse a phip xarray dataset by some group in the metadata.
    """

    # Check to see if the group(s) is/are available
    groups_avail = ds[f"{collapse_dim}_metadata"].values
    for group in by:
        if group not in groups_avail:
            raise KeyError(
                f"{group} is not included as a column in the {collapse_dim} table. The available groups are {groups_avail}"
            )

    # grab relavent annotation tables
    collapse_df = ds[f"{collapse_dim}_table"].to_pandas()

    # Create group-able dataset by assigning table columns to a coordinate
    if len(by) == 1:
        coord = collapse_df[by[0]]
        coord_ds = ds.assign_coords(coord=(f"{collapse_dim}_id", coord))

    # if were grouping by multiple things, we need to zip 'em into a tuple coord
    else:
        common_dim = f"{collapse_dim}_id"
        coor_arr = np.empty(len(ds[common_dim]), dtype=object)
        coor_arr[:] = list(zip(*(collapse_df[f].values for f in by)))
        coord_ds = ds.assign_coords(
            coord=xr.DataArray(coor_arr, collapse_dims=common_dim)
        )

    # Save dat memory, also, we will perform custom grouping's on the annotation tables
    del coord_ds["sample_table"]
    del coord_ds["peptide_table"]

    # Perform the reduction on all data tables.
    collapsed_enrichments = coord_ds.groupby("coord").reduce(agg_func).transpose()

    # Once the data tables are grouped we have no use for first copy.
    del coord_ds

    # Compile each of the collapsed xarray variables.
    collapsed_xr_dfs = {
        f"{dt}": (
            ["peptide_id", "sample_id"],
            collapsed_enrichments[f"{dt}"].to_pandas(),
        )
        for dt in set(list(collapsed_enrichments.data_vars))
    }

    # get collapsed table
    # TODO, you could also offer a "first()" option here.
    collapsed_anno_table = throw_mismatched_features(collapse_df, by)

    # Compute mean pairwise correlation for all groups,
    # for all enrichment layers - and add it to the
    # resulting collapsed sample table
    # TODO, this could be slightly expensive, and possibly more
    # well rounded.
    # NOTE TODO
    # if compute_pw_cc:

    dims = set(["sample", "peptide"])
    fixed_dim = list(dims - set([collapse_dim]))[0]

    # Insert the correct annotation tables to out dictionary of variables
    # NOTE AnnoTable.
    collapsed_xr_dfs[f"{collapse_dim}_table"] = (
        [f"{collapse_dim}_id", "{dim}_metadata"],
        collapsed_anno_table,
    )

    collapsed_xr_dfs[f"{fixed_dim}_table"] = (
        ["{fixed_dim}_id", "{fixed_dim}_metadata"],
        ds[f"{fixed_dim}_table"].to_pandas(),
    )

    pds = xr.Dataset(
        collapsed_xr_dfs,
        coords={
            "sample_id": collapsed_xr_dfs["sample_table"][1].index.values,
            "peptide_id": collapsed_xr_dfs["peptide_table"][1].index.values,
            "sample_metadata": collapsed_xr_dfs["sample_table"][1].columns.values,
            "peptide_metadata": collapsed_xr_dfs["peptide_table"][1].columns.values,
        },
    )
    pds.attrs["collapsed_group"] = group
    return pds
