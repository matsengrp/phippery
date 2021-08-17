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
    keep only features which are shared by

    This is similar to
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


def collapse_groups(ds, by, table="sample", agg_func=np.mean, compute_pw_cc=False):
    """
    Collapse a phip xarray dataset by some group in the metadata.
    """

    # Check to see if the group(s) is/are available
    groups_avail = ds[f"{table}_metadata"].values
    for group in by:
        if group not in groups_avail:
            raise KeyError(
                f"{group} is not included as a column in the {table} table. The available groups are {groups_avail}"
            )

    # grab relavent annotation tables
    collapse_df = ds[f"{table}_table"].to_pandas()
    tables = set(["sample_table", "peptide_table"])
    fixed_var = list(tables - set([f"{table}_table"]))[0]

    # Create group-able dataset by assigning table columns to a coordinate
    if len(by) == 1:
        coord = collapse_df[by[0]]
        coord_ds = ds.assign_coords(coord=(f"{table}_id", coord))

    # if were grouping by multiple things, we need to zip 'em into a tuple coord
    else:
        common_dim = f"{table}_id"
        coor_arr = np.empty(len(ds[common_dim]), dtype=object)
        coor_arr[:] = list(zip(*(collapse_df[f].values for f in by)))
        coord_ds = ds.assign_coords(coord=xr.DataArray(coor_arr, dims=common_dim))

    # save dat memory
    del coord_ds["sample_table"]
    del coord_ds["peptide_table"]
    print(coord_ds)
    print("YOOOOOOOOOOOOOOOOOOOO")

    # perform the reduction
    collapsed_enrichments = coord_ds.groupby("coord").reduce(agg_func).transpose()

    # oce grouped, we have no use for first copy, delete
    del coord_ds

    # compile each of the collapsed xarray variables
    collapsed_xr_dfs = {
        f"{dt}": (
            ["peptide_id", "sample_id"],
            collapsed_enrichments[f"{dt}"].to_pandas(),
        )
        for dt in set(list(collapsed_enrichments.data_vars))
    }

    # get collapsed table
    # TODO, you could also offer a "First" option here.
    csm = throw_mismatched_features(collapse_df, by)

    # Compute mean pairwise correlation for all groups,
    # for all enrichment layers - and add it to the
    # resulting collapsed sample table
    # TODO, this could be
    if compute_pw_cc:
        for data_table in list(
            set(ds.data_vars) - set(["sample_table", "peptide_table"])
        ):
            pw_cc = pairwise_correlation_by_sample_group(ds, by, data_table)
            assert len(pw_cc) == len(csm)
            csm = csm.merge(pw_cc, left_index=True, right_index=True)

    # insert the correct annotation tables to
    # NOTE AnnoTable.
    collapsed_xr_dfs[f"{table}_table"] = (
        [f"{table}_id", "sample_metadata"],
        csm,
    )

    collapsed_xr_dfs[f"{fixed_var}"] = (
        ["peptide_id", "peptide_metadata"],
        ds[f"{fixed_var}"].to_pandas(),
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


def pairwise_correlation_by_group(
    ds, group="sample_ID", data_table="counts", column_prefix=None
):
    """
    a method which computes pairwise cc for all
    sample in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of samples
    in the group.
    """

    if group not in ds.sample_metadata.values:
        raise ValueError("{group} does not exist in sample metadata")

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


####################################################
####################################################


def collapse_sample_groups(ds, group, agg_func=np.mean, compute_pw_cc=False):
    """
    Collapse a phip xarray dataset by some group in the metadata.
    """

    if group not in ds.sample_metadata.values:
        raise KeyError(
            f"{group} is not included as a column in the sample table. The available groups are {ds.sample_metadata.values}"
        )

    try:
        coord = ds.sample_table.loc[:, group].values.astype(int)
    except ValueError:
        raise ValueError(
            f"All factor level values in {group} must be able to convert to int"
        )

    coord_ds = ds.assign_coords(coord=("sample_id", coord))
    del coord_ds["sample_table"]
    del coord_ds["peptide_table"]
    collapsed_enrichments = coord_ds.groupby("coord").reduce(agg_func).transpose()

    collapsed_xr_dfs = {
        f"{dt}": (
            ["peptide_id", "sample_id"],
            collapsed_enrichments[f"{dt}"].to_pandas(),
        )
        for dt in set(list(collapsed_enrichments.data_vars))
    }

    sample_table_df = ds.sample_table.to_pandas()
    collapsed_sample_metadata = defaultdict(list)
    for i, (tech_rep_id, tech_rep_meta) in enumerate(sample_table_df.groupby(group)):
        for column, value in tech_rep_meta.iteritems():
            v = value.values
            if np.all(v == v[0]) or np.all([n != n for n in v]):
                collapsed_sample_metadata[column].append(v[0])

    to_throw = [
        key for key, value in collapsed_sample_metadata.items() if len(value) < i + 1
    ]

    [collapsed_sample_metadata.pop(key) for key in to_throw]

    csm = pd.DataFrame(collapsed_sample_metadata)
    csm[group] = csm[group].astype(int)
    csm.set_index(group, inplace=True)

    if compute_pw_cc:
        for data_table in list(
            set(ds.data_vars) - set(["sample_table", "peptide_table"])
        ):
            pw_cc = pairwise_correlation_by_sample_group(ds, group, data_table)
            assert len(pw_cc) == len(csm)
            csm = csm.merge(pw_cc, left_index=True, right_index=True)

    collapsed_xr_dfs["sample_table"] = (
        ["sample_id", "sample_metadata"],
        csm,
    )

    collapsed_xr_dfs["peptide_table"] = (
        ["peptide_id", "peptide_metadata"],
        ds.peptide_table.to_pandas(),
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


def pairwise_correlation_by_sample_group(
    ds, group="sample_ID", data_table="counts", column_prefix=None
):
    """
    a method which computes pairwise cc for all
    sample in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of samples
    in the group.
    """

    if group not in ds.sample_metadata.values:
        raise ValueError("{group} does not exist in sample metadata")

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
