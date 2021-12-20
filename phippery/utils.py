"""
@File: utils.py

@Author: Jared Galloway

This file will include some helpful functions
for the phippery package CLI. The primary
data struct we put counts data into is
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
from collections import defaultdict

from phippery.phipdata import get_sample_table
from phippery.phipdata import get_peptide_table
from phippery.phipdata import get_annotation_table


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


def iter_groups(ds, by, dim="sample"):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by items in the metadata of either dimension.
    """

    table = copy.deepcopy(ds[f"{dim}_table"].to_pandas().convert_dtypes())
    for group, group_df in table.groupby(by):
        group_ds = ds.loc[{f"{dim}_id": list(group_df.index.values)}]
        # group_ds = ds.loc[dict(peptide_id=list(group_df.index.values))]
        yield group, group_ds


def iter_sample_groups(*args):
    """
    DEPRECATED - use 'iter_groups()' instead

    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the sample metadata coodinate.
    """

    # sample_table = ds.sample_table.to_pandas()
    # for group, group_st in sample_table.groupby(groupby):
    #    group_ds = ds.loc[dict(sample_id=list(group_st["sample_id"].values))]
    #    yield group, group_ds

    return iter_groups(*args, "sample")


def iter_peptide_groups(*args):

    """DEPRECATED - use 'iter_groups()' instead

    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the peptide metadata coodinate. """

    # peptide_table = ds.peptide_table.to_pandas()
    # for group, group_st in peptide_table.groupby(groupby):
    #    group_ds = ds.loc[dict(peptide_id=list(group_st.index.values))]
    #    yield group, group_ds

    return iter_groups(*args, "peptide")


def id_coordinate_from_query(ds, query_df):

    """

    """

    # TODO make this 'dim', instead of 'Type'
    # TODO run checks, Raise Errors

    # st = ds.sample_table.to_pandas().infer_objects()
    sq = list(query_df.loc[query_df["dimension"] == "sample", "expression"].values)
    sid = sample_id_coordinate_from_query(ds, sq)

    # pt = ds.peptide_table.to_pandas().infer_objects()
    pq = list(query_df.loc[query_df["dimension"] == "peptide", "expression"].values)
    pid = peptide_id_coordinate_from_query(ds, pq)

    return sid, pid


def peptide_id_coordinate_from_query(ds, query_list: list, *args, **kwargs):
    """Take in a list of queries and return the peptide id index resulting
    from query """

    if len(query_list) == 0:
        return list(ds.peptide_id.values)

    peptide_table = get_peptide_table(ds)
    return list(peptide_table.query(" & ".join(query_list)).index.values)


def sample_id_coordinate_from_query(ds, query_list: list, *args, **kwargs):
    """Take in a list of queries and return the sample id index resulting
    from query """

    if len(query_list) == 0:
        return list(ds.sample_id.values)

    sample_table = get_sample_table(ds)
    return list(sample_table.query(" & ".join(query_list)).index.values)
