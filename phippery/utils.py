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


def sample_id_coordinate_subset(
    ds,
    where,
    is_equal_to=None,
    is_not_equal_to=None,
    is_greater_than=None,
    is_less_than=None,
    is_in=None,
    is_valid=None,
):
    """
    grab the subset of sample id's representing a factor
    group in the sample metadata

    Only one of the keyword arguments can be provided.

    :param: "is_equal_to" expects a single value, giving all sample_id's
        where "where" is equal to the given value.

    :param: ""
    """

    # TODO Make this a general function
    # for sample and peptide table
    # you're going to need to figure out the
    table = "sample_table"

    num_kw_args = [
        bool(arg)
        for arg in [is_equal_to, is_greater_than, is_less_than, is_in, is_valid]
    ]

    if sum(num_kw_args) > 1:
        raise ValueError("You may only provide one of the keyword arguments")

    if is_equal_to:
        return list(
            ds.sample_id.where(
                ds[table].loc[:, where] == is_equal_to, drop=True
            ).sample_id.values
        )

    elif is_not_equal_to:
        return list(
            ds.sample_id.where(
                ds[table].loc[:, where] != is_not_equal_to, drop=True
            ).sample_id.values
        )

    elif is_greater_than:
        return list(
            ds.sample_id.where(
                ds[table].loc[:, where] >= is_greater_than, drop=True
            ).sample_id.values
        )

    elif is_less_than:
        return list(
            ds.sample_id.where(
                ds[table].loc[:, where] < is_less_than, drop=True
            ).sample_id.values
        )

    elif is_in:
        return list(
            ds.sample_id.where(
                ds[table].loc[:, where].isin(is_in), drop=True
            ).sample_id.values
        )

    else:
        return list(
            ds[table]
            .where(ds[table].loc[:, where] == ds[table].loc[:, where], drop=True,)
            .sample_id.values
        )


def peptide_id_coordinate_subset(
    ds,
    where,
    is_equal_to=None,
    is_not_equal_to=None,
    is_greater_than=None,
    is_less_than=None,
    is_in=None,
    is_valid=None,
):
    """
    grab the subset of peptide id's representing a factor
    group in the peptide metadata

    Only one of the keyword arguments can be provided.
    """

    # TODO Make this a general function
    # for peptide and peptide table
    # you're going to need to figure out the
    table = "peptide_table"

    num_kw_args = [
        bool(arg)
        for arg in [is_equal_to, is_greater_than, is_less_than, is_in, is_valid]
    ]

    if sum(num_kw_args) > 1:
        raise ValueError("You may only provide one of the keyword arguments")

    if is_equal_to:
        return list(
            ds.peptide_id.where(
                ds[table].loc[:, where] == is_equal_to, drop=True
            ).peptide_id.values
        )

    elif is_not_equal_to:
        return list(
            ds.peptide_id.where(
                ds[table].loc[:, where] != is_not_equal_to, drop=True
            ).peptide_id.values
        )

    elif is_greater_than:
        return list(
            ds.peptide_id.where(
                ds[table].loc[:, where] >= is_greater_than, drop=True
            ).peptide_id.values
        )

    elif is_less_than:
        return list(
            ds.peptide_id.where(
                ds[table].loc[:, where] < is_less_than, drop=True
            ).peptide_id.values
        )

    elif is_in:
        return list(
            ds.peptide_id.where(
                ds[table].loc[:, where].isin(is_in), drop=True
            ).peptide_id.values
        )

    else:
        return list(
            ds[table]
            .where(ds[table].loc[:, where] == ds[table].loc[:, where], drop=True,)
            .peptide_id.values
        )


def sample_id_subset(ds, where, is_equal_to):
    """
    grab the subset of sample id's representing a factor
    group in the peptide metadata
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
