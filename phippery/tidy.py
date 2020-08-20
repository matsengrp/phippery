"""
@File: utils.py

@Author: Jared Galloway

This file will include some helpful functions
for the phippery package CLI. The primary
data struct we put counts data into is

So far it includes functions to"

* compile tsv files into a phip dataset
* TODO check phip_dataset attributed

"""

# local
from phippery.utils import iter_sample_groups

# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as st

# built-in python3
import os
import copy
import itertools
from functools import reduce


def tidy_ds(ds):
    """
    return the provided dataset in 'tall form'.

    This means that all information from the entire dataset and
    all the tables it includes will be 'melted' into a single dataframe.
    This format is far less effecient in terms of storage, and should not
    be used on the entire dataset, but rather, some subset of the dataset
    when you are ready to plot the relevent information.

    The means that for each row, you will get the following columns:

    sample_id,
    peptide_id,
    all sample metadata (a column for each),
    all peptide metadata (a column for each),
    all enrichment tables (a column for each).

    Ideal for gg plotting.
    """

    # melt all data tables in the dataset
    melted_data = [
        ds[f"{dt}"]
        .to_pandas()
        .reset_index()
        .melt(id_vars=["peptide_id"])
        .rename({"value": f"{dt}"}, axis=1)
        for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    ]

    # merge all melted tables so each get's a column in the final df
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, on=["peptide_id", "sample_id"]), melted_data
    )

    # grab only the columns of the metadata we want in the resulting dataframe
    peptide_table = ds.peptide_table.to_pandas().reset_index()
    sample_table = ds.sample_table.to_pandas().reset_index()

    # merge the metadata into the melted datatables
    data_peptide = merged_counts_df.merge(peptide_table, on="peptide_id")
    return data_peptide.merge(sample_table, on="sample_id")


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
    # for group, group_ds in ds.groupby(ds.sample_table.loc[:, group]):
    for group, group_ds in iter_sample_groups(ds, group):
        groups.append(group)
        n.append(len(group_ds.sample_id.values))
        if len(group_ds.sample_id.values) < 2:
            pw_cc.append(0)
            continue
        correlations = []
        for sample_ids in itertools.combinations(group_ds.sample_id.values, 2):
            sample_0_enrichment = data.loc[:, sample_ids[0]]
            sample_1_enrichment = data.loc[:, sample_ids[1]]
            correlation = (
                st.pearsonr(sample_0_enrichment, sample_1_enrichment)[0]
                if np.any(sample_0_enrichment != sample_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(round(sum(correlations) / len(correlations), 3))

    return pd.DataFrame({"group": groups, f"{data_table}_pw_cc": pw_cc, "n_reps": n})


def tidy_sample_table(ds):
    """simply return the sample df"""
    return copy.deepcopy(ds.sample_table.to_pandas())


def tidy_peptide_table(ds):
    """simply return the sample df"""
    return copy.deepcopy(ds.peptide_table.to_pandas())
