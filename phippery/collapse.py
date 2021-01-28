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


def collapse_sample_groups(
    ds, group, agg_func=lambda x: np.mean(x, axis=1), compute_pw_cc=False
):
    """
    Collapse a phip xarray dataset by some group in the metadata.
    """

    if group not in ds.sample_metadata.values:
        raise KeyError(
            f"{group} is not included as a column in the sample table. The available groups are {ds.sample_metadata.values}"
        )

    coord = ds.sample_table.loc[:, group].values.astype(int)
    coord_ds = ds.assign_coords(coord=("sample_id", coord))
    collapsed_enrichments = (
        coord_ds.groupby("coord").map(lambda x: agg_func(x),).transpose()
    )

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

    csm = pd.DataFrame(collapsed_sample_metadata).set_index(group)
    if compute_pw_cc:
        pw_cc = pairwise_correlation_by_sample_group(ds, group)
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


def pairwise_correlation_by_sample_group(ds, group="sample_ID", data_table="counts"):
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
        groups.append(s_group)
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

        ###############
        # if len(group_ds["sample_id"].values) != 2:
        #    pw_cc.append(1.0)
        #    continue
        # sample_0_enrichment = data.loc[:, group_ds["sample_id"].values[0]]
        # sample_1_enrichment = data.loc[:, group_ds["sample_id"].values[1]]

        # correlation = (
        #    st.pearsonr(sample_0_enrichment, sample_1_enrichment)[0]
        #    if np.any(sample_0_enrichment != sample_1_enrichment)
        #    else 1.0
        # )
        # pw_cc.append(round(correlation, 5))
        ################

    ret = pd.DataFrame(
        {f"{group}": groups, f"{group}_pw_cc": pw_cc, f"{group}_n_reps": n}
    ).set_index(group)

    return ret
