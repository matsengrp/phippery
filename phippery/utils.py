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


def trim_index(sequence):
    return "".join([nt for nt in sequence if nt.isupper()])


def convert_peptide_metadata_to_fasta(peptide_metadata, out):
    """
    Take in peptide metadata dataframe, and write a fasta
    format representation of the oligos
    """

    fasta_fp = open(out, "w")
    peptide_metadata = pd.read_csv(peptide_metadata, index_col=0, header=0)
    requirements = ["Oligo"]
    assert peptide_metadata.index.name == "peptide_id"
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


# This could be generalized to do some helpful things with the enirchment tables
def id_coordinate_subset(
    ds,
    where,
    table="sample_table",
    is_equal_to=None,
    is_not_equal_to=None,
    is_greater_than=None,
    is_greater_than_or_equal_to=None,
    is_less_than=None,
    is_less_than_or_equal_to=None,
    is_in=None,
    is_valid=None,
):
    """
    a general function to compute the coordinate dimensions given some conditions.
    """

    if table not in ["sample_table", "peptide_table"]:
        raise ValueError(
            f"{table} is not a valid data table for {ds}\n Available data tables are: 'sample_table' or 'peptide_table'"
        )

    if table == "sample_table":
        metadata = "sample_metadata"
        metadata_features = ds[metadata]
        coord = "sample_id"
    else:
        metadata = "peptide_metadata"
        metadata_features = ds[metadata]
        coord = "peptide_id"

    if where not in metadata_features:
        raise ValueError(
            f"{where} is not in the sample metadata\n Available options are: {metadata_features.values}"
        )

    num_kw_args = [
        0 if arg is None else 1
        for arg in [
            is_equal_to,
            is_not_equal_to,
            is_greater_than,
            is_greater_than_or_equal_to,
            is_less_than,
            is_less_than_or_equal_to,
            is_in,
            is_valid,
        ]
    ]

    if sum(num_kw_args) != 1:
        raise ValueError(
            "You must provide exactly one of the keyword conditional arguments"
        )

    table = ds[table]
    dim = table.loc[{metadata: where}]
    coordinate_ids = ds[coord]

    if is_equal_to is not None:
        return coordinate_ids[dim == is_equal_to].values

    elif is_not_equal_to is not None:
        return coordinate_ids[dim != is_not_equal_to].values

    elif is_greater_than is not None:
        return coordinate_ids[dim > is_greater_than].values

    elif is_greater_than_or_equal_to is not None:
        return coordinate_ids[dim >= is_greater_than_or_equal_to].values

    elif is_less_than is not None:
        return coordinate_ids[dim < is_less_than].values

    elif is_less_than_or_equal_to is not None:
        return coordinate_ids[dim <= is_less_than_or_equal_to].values

    elif is_in is not None:
        return coordinate_ids[dim.isin(is_in)].values

    else:
        return coordinate_ids[dim == dim].values


def sample_id_coordinate_subset(ds, where, **kwargs):
    return id_coordinate_subset(ds, where, "sample_table", **kwargs)


def peptide_id_coordinate_subset(ds, where, **kwargs):
    return id_coordinate_subset(ds, where, "peptide_table", **kwargs)
