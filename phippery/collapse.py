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

# built-in python3
from collections import defaultdict


def collapse_sample_groups(ds, group, agg_func=lambda x: np.mean(x, axis=1)):
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

    collapsed_xr_dfs["sample_table"] = (
        ["sample_id", "sample_metadata"],
        pd.DataFrame(collapsed_sample_metadata).set_index(group),
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
    # pds.attrs["sample_coord_dim"] = group
    # pds.attrs["peptide_coord_dim"] = "peptide_id"
    return pds
