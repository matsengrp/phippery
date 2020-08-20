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


def convert_peptide_metadata_to_fasta(peptide_metadata, out):
    """
    Take in peptide metadata dataframe, and write a fasta
    format representation of the oligos
    """

    def trim_index(sequence):
        return "".join([nt for nt in sequence if nt.isupper()])

    fasta_fp = open(out, "w")
    peptide_metadata = pd.read_csv(peptide_metadata, index_col=0, header=0)
    requirements = ["Oligo"]
    assert peptide_metadata.index.name == "ID"
    assert np.all([x in peptide_metadata.columns for x in requirements])
    for index, row in peptide_metadata.iterrows():
        ref_sequence = trim_index(row["Oligo"])
        fasta_fp.write(f">{index}\n{ref_sequence}\n")


def get_all_sample_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a sample table column
    """

    all_exp = ds.sample_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def get_all_peptide_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a peptide table column
    """

    all_exp = ds.peptide_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def melted_mean_collapse_groups(ds, by="sample_ID"):
    """
    This function melts a dataset after taking the average of
    all values in sample metadata coordinate group.
    The nature of this function forces us to throw out sample
    table information for all but the first of the group
    encountered in the dataset - defined by the "by" category.
    """

    # TODO This is ugly and should re-structure the
    # the entire dataset? Or maybe not?
    # Either way, even if this continued as a melt/mean method,
    # let's clean it up

    sample_table = ds.sample_table.to_pandas()
    # create a "by" sample table, which only takes the first instance
    # of a group in the sample metadata, and sets that as the index.

    # grab the first sample in each group to keep it's sample metadata.
    # This metadata will be representative of all samples which get averaged.
    to_keep = [
        group.sample_id.values[0] for group_f, group in iter_sample_groups(ds, by)
    ]
    bio_sample_table = sample_table.loc[to_keep, :].set_index(by, verify_integrity=True)

    # a new counts table, where each of the groups are averaged, then the table is melted
    avg_bio_reps_df = ds.counts.groupby(ds.sample_table.loc[:, by]).mean().to_pandas()
    melted = avg_bio_reps_df.reset_index().melt(id_vars=["peptide_id"])
    melted.rename({"sample_table": by}, axis=1, inplace=True)

    # just our peptide table - it's perfect how it is except we need the pid to merge by.
    peptide_table = ds.peptide_table.to_pandas().reset_index()

    melted_peptide = melted.merge(peptide_table, on="peptide_id")
    all_melted = melted_peptide.merge(bio_sample_table, on=by)

    return all_melted


def melt_ds(ds, sample_metadata_columns="all", peptide_metadata_columns="all"):
    """
    return an new dataset, with all data compltely melted.
    This means that instead of a counts matrix
    Ideal for gg plotting.
    """

    melted_data = [
        ds[f"{dt}"]
        .to_pandas()
        .reset_index()
        .melt(id_vars=["peptide_id"])
        .rename({"value": f"{dt}"}, axis=1)
        for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    ]
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, on=["peptide_id", "sample_id"]), melted_data
    )
    peptide_table = ds.peptide_table.to_pandas().reset_index()
    sample_table = ds.sample_table.to_pandas().reset_index()
    if sample_metadata_columns != "all":
        sample_table = sample_table[["sample_id"] + sample_metadata_columns]
    if peptide_metadata_columns != "all":
        peptide_table = peptide_table[["peptide_id"] + peptide_metadata_columns]
    data_peptide = merged_counts_df.merge(peptide_table, on="peptide_id")
    return data_peptide.merge(sample_table, on="sample_id")


def iter_sample_groups(ds, column):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the sample metadata coodinate.
    """
    for group, group_ds_idx in ds.sample_id.groupby(ds.sample_table.loc[:, column]):
        group_ds = ds.loc[dict(sample_id=group_ds_idx.sample_id.values)]
        yield group, group_ds


def iter_peptide_groups(ds, column):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the peptide metadata coodinate.
    """
    for group, group_ds_idx in ds.peptide_id.groupby(ds.peptide_table.loc[:, column]):
        group_ds = ds.loc[dict(peptide_id=group_ds_idx.peptide_id.values)]
        yield group, group_ds


def sample_id_subset(ds, where, is_equal_to):
    """
    grab the subset of sample id's representing a factor
    group in the sample metadata
    """
    return list(
        ds.sample_id.where(
            ds.sample_table.loc[:, where] == is_equal_to, drop=True
        ).sample_id.values
    )


def peptide_id_subset(ds, where, is_equal_to):
    """
    grab the subset of sample id's representing a factor
    group in the peptide metadata
    """
    return list(
        ds.peptide_id.where(
            ds.peptide_table.loc[:, where] == is_equal_to, drop=True
        ).peptide_id.values
    )


def valid_factor_peptide_id_subset(ds, where_is_valid):
    """
    return the subset of peptide id's which are non NaN
    """

    return list(
        ds.peptide_id.where(
            ds.peptide_table.loc[:, where_is_valid]
            == ds.peptide_table.loc[:, where_is_valid],
            drop=True,
        ).peptide_id.values
    )


def valid_factor_sample_id_subset(ds, where_is_valid):
    """
    return the subset of peptide id's which are non NaN
    """

    return list(
        ds.sample_id.where(
            ds.sample_table.loc[:, where_is_valid]
            == ds.sample_table.loc[:, where_is_valid],
            drop=True,
        ).sample_id.values
    )


def subset_ds_by_sample_factor(ds, where, is_equal_to):
    """
    grab the subset loc of the ds representing a factor
    group in the sample metadata
    """
    sample_ids = list(
        ds.sample_id.where(
            ds.sample_table.loc[:, where] == is_equal_to, drop=True
        ).sample_id.values
    )

    return ds.loc[dict(sample_id=sample_ids)]


def subset_ds_by_peptide_factor(ds, where, is_equal_to):
    """
    grab the subset loc of the ds representing a factor
    group in the peptide metadata
    """
    peptide_ids = list(
        ds.peptide_id.where(
            ds.peptide_table.loc[:, where] == is_equal_to, drop=True
        ).peptide_id.values
    )

    return ds.loc[dict(peptide_id=peptide_ids)]
