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
from functools import reduce
import itertools

# local
from phippery.utils import iter_sample_groups
from phippery.utils import iter_groups


# TODO this needs to be tested
def throw_mismatched_features(df, by):

    """When you collapse by some set of columns in the dataframe,
    keep only features which homogeneous within groups.

    This is similar to 'DataFrameGroupby.first()' method,
    but instead of keeping the first factor level appearing for each group
    category, we only throw any features which are note homogeneous within groups.
    You could achieve the same functionality by listing features you know to be
    homogeneous in the 'by' parameter."""

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


def mean_pw_cc_by_multiple_tables(ds, by, dim="sample", data_tables="all"):

    """A wrapper for computing pw cc within groups defined
    with the 'by' parameter. Merges the correlations into a
    single table"""

    # Compute mean pw cc on all possible data tables
    if data_tables == "all":
        data_tables = list(set(ds.data_vars) - set(["sample_table", "peptide_table"]))

    # Some error handling
    if dim not in ["sample", "peptide"]:
        raise ValueError(f"parameter 'dim' must be either 'sample' or 'peptide'")

    groups_avail = ds[f"{dim}_metadata"].values
    for data_table in data_tables:
        if data_table not in ds:
            raise KeyError(f"{data_table} is not included in dataset.")

    for group in by:
        if group not in groups_avail:
            raise KeyError(
                f"{group} is not included as a column in the {dim} table. The available groups are {groups_avail}"
            )

    # Compute mean pw cc on all data layers - resulting in a df for each
    corr_dfs = [
        mean_pw_cc_by(ds, by, data_table=data_table, dim=dim)
        for data_table in data_tables
    ]

    # return a single merged df containing info for all data layer pw cc
    return reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        corr_dfs,
    )


# TODO - this needsd to be generalized as well.
def mean_pw_cc_by(ds, by, data_table="counts", dim="sample"):

    """Computes pairwise cc for all
    dim in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of dims
    in the group."""

    data = ds[f"{data_table}"].to_pandas()

    # TODO Let's allocate mem right here instead of hefty appending.
    groups, pw_cc, n = [], [], []

    for s_group, group_ds in iter_groups(ds, by, dim):
        groups.append(int(s_group))
        n.append(len(group_ds[f"{dim}_id"].values))

        if len(group_ds[f"{dim}_id"].values) < 2:
            pw_cc.append(1.0)
            continue

        # TODO Same as above
        correlations = []
        for dim_ids in itertools.combinations(group_ds[f"{dim}_id"].values, 2):
            dim_0_enrichment = data.loc[:, dim_ids[0]]
            dim_1_enrichment = data.loc[:, dim_ids[1]]
            # TODO do I need to do this?
            correlation = (
                st.pearsonr(dim_0_enrichment, dim_1_enrichment)[0]
                if np.any(dim_0_enrichment != dim_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(round(sum(correlations) / len(correlations), 5))

    name = "_".join(by)
    column_prefix = f"{name}_{data_table}"

    ret = pd.DataFrame(
        {
            f"{name}": groups,
            f"{column_prefix}_pw_cc": np.array(pw_cc).astype(np.float64),
            f"{column_prefix}_n_reps": np.array(n).astype(np.int),
        }
    ).set_index(name)

    return ret


def collapse_sample_groups(*args, **kwargs):
    """
    DEPRECATED - SEE COLLAPSE GROUPS
    """
    return collapse_groups(*args, **kwargs, collapse_dim="sample")


def collapse_groups(
    ds, by, collapse_dim="sample", agg_func=np.mean, compute_pw_cc=False, **kwargs
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

    # define collapse and fixed df
    dims = set(["sample", "peptide"])
    fixed_dim = list(dims - set([collapse_dim]))[0]

    # grab relavent annotation tables
    # TODO infer objects, here, right?
    # Question is, is there any way that groupby changes when you have
    # different datatypes when compared to "objects"
    # TODO write a unit test for this ^^.
    collapse_df = ds[f"{collapse_dim}_table"].to_pandas()
    fixed_df = ds[f"{fixed_dim}_table"].to_pandas()

    # Create group-able dataset by assigning table columns to a coordinate
    # TODO How do we re-name this according to the thing we're collapsing on???
    if len(by) == 1:
        coord = collapse_df[by[0]]
        coord_ds = ds.assign_coords({f"{by[0]}": (f"{collapse_dim}_id", coord)})
    else:
        print(f"WARNING: Nothing available, here")
        return None

    # TODO
    # if were grouping by multiple things, we need to zip 'em into a tuple coord
    # psuedo-code
    # else:
    #    common_dim = f"{collapse_dim}_id"
    #    coor_arr = np.empty(len(ds[common_dim]), dtype=object)
    #    coor_arr[:] = list(zip(*(collapse_df[f].values for f in by)))
    #    coord_ds = ds.assign_coords(
    #        coord=xr.DataArray(coor_arr, collapse_dims=common_dim)
    #    )

    # Save dat memory, also, we will perform custom grouping's on the annotation tables
    # so, no need for them here
    del coord_ds["sample_table"]
    del coord_ds["peptide_table"]
    del coord_ds["sample_metadata"]
    del coord_ds["peptide_metadata"]

    # Perform the reduction on all data tables.
    collapsed_enrichments = coord_ds.groupby(f"{by[0]}", squeeze=False).reduce(agg_func)

    if collapse_dim == "sample":
        collapsed_enrichments = collapsed_enrichments.transpose()

    # Once the data tables are grouped we have no use for first copy.
    # TODO can we do this in place?
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
    # TODO, you could also offer a "first" option here.
    # Collapsed Annotation Table, 'cat' for brevity
    cat = throw_mismatched_features(collapse_df, by)
    # print(cat)

    # Compute mean pairwise correlation for all groups,
    # for all enrichment layers - and add it to the
    # resulting collapsed sample table
    # TODO,
    if compute_pw_cc:
        print(f"WARNING - ONLY ON SAMPLES ATM")
        mean_pw_cc = mean_pw_cc_by(ds, by, **kwargs)
        cat = cat.merge(mean_pw_cc, left_index=True, right_index=True)

    # Insert the correct annotation tables to out dictionary of variables
    # NOTE AnnoTable.
    collapsed_xr_dfs[f"{collapse_dim}_table"] = (
        [f"{collapse_dim}_id", f"{collapse_dim}_metadata"],
        cat,
    )

    collapsed_xr_dfs[f"{fixed_dim}_table"] = (
        [f"{fixed_dim}_id", f"{fixed_dim}_metadata"],
        fixed_df,
    )

    pds = xr.Dataset(
        collapsed_xr_dfs,
        coords={
            f"{collapse_dim}_id": collapsed_xr_dfs[f"{collapse_dim}_table"][
                1
            ].index.values,
            f"{fixed_dim}_id": collapsed_xr_dfs[f"{fixed_dim}_table"][1].index.values,
            f"{collapse_dim}_metadata": collapsed_xr_dfs[f"{collapse_dim}_table"][
                1
            ].columns.values,
            f"{fixed_dim}_metadata": collapsed_xr_dfs[f"{fixed_dim}_table"][
                1
            ].columns.values,
        },
    )
    return pds
