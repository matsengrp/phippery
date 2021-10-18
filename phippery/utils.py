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


#def string_feature(ds, feature: str, verbosity: str):
#    """
#    Take in a single feature
#    """
#    ser = st.session_state.sample_table[sample_sum]
#    dt = ser.dtype
#    des = ser.describe()
#    if dt == pd.StringDtype():
#
#        # think about this
#        levels = [f for f in list(set(ser)) if f == f]
#        num_levels = len(levels)
#        if num_levels >= 50:
#            st.write(f"""
#            {sample_sum} is a String feature with too many factor levels (n = {num_levels}) to describe, 
#            here (probably a desciptive string, and thus, not great visualization feature)
#
#            Summary:
#
#            {des}
#            """)
#        elif num_levels < 50 and num_levels >=2:
#            all_factors = ", ".join(levels)
#            st.write(f"""
#            {sample_sum} is a categorical feature with the following factor levels to investigate, 
#
#            **{all_factors}**
#
#            Some example query statements:
#
#            > {sample_sum} in ['{levels[0]}', '{levels[1]}', ...]
#
#            > {sample_sum} not in ['{levels[0]}', '{levels[-1]}', ...]
#
#            > {sample_sum} != '{levels[-2]}'
#
#            """)
#
#            f"""
#            ### {sample_sum} Factor level counts
#            """
#            st.write(ser.value_counts())
#
#            f"""
#            ### {sample_sum} Factor level counts
#            """
#            st.write(des)
#        else:
#            st.write(f"""
#                There's only a single factor level, {levels[0]}, across all samples.
#            """
#            )
#
#    elif dt == pd.BooleanDtype():
#        st.write(f"""
#            Boolean Feature:
#
#            Some example query statements:
#
#            > {sample_sum} == True
#
#            > {sample_sum} == False
#        """
#        )
#        f"""
#        ### {sample_sum} Factor level counts
#        """
#        st.write(ser.value_counts())
#    elif dt == pd.Int64Dtype() or dt == pd.Float64Dtype():
#        st.write(f"""
#            Numerical Feature:
#
#            Some example query statements:
#
#            > {sample_sum} >= {des[1]}
#
#        """
#        )
#        f"""
#        ### {sample_sum} Summary
#        """
#        des = ser.describe()
#        st.write(des)
#    else:
#        st.error("Never seen this Dtype before!")
#
#
#def string_ds(ds, verbosity: int):
#    """
#    Summarize the data in a given dataset
#
#    If verbosity flag is set to zero, this will print the 
#    basic information about number of
#    samples, peptides, and enrichment layers
#    in a given dataset. With a verbosity of one (-v) you will
#    get also information about annotations and available datatypes.
#    If verbosity flag is set to two (-vv) - Print 
#    detailed information about all data tables including annotation
#    feature types and distribution of enrichment values for each layer.
#    A verbosity of three will basically loop through all features
#    """
#
#    # initialize formatting for each of the three major facets of a dataset
#    table_strings = {
#        'sample_table' : """
#Sample Table:
#-------------
#""",
#
#        'peptide_table' : """
#Peptide Table:
#--------------
#""",
#
#        'enrichments' : """
#Enrichment Matrices:
#--------------------
#"""
#    }
#
#    for dimension in ["sample", "peptide"]:
#
#        num_dimensions= len(ds[f"{dimension}_id"].values)
#        dimension_annotations = list(ds[f"{dimension}_metadata"].values)
#        
#        dimension_annotation_strings = {
#            f'{l}':f"""
#            * {l}
#            """
#        for l in dimension_annotations
#        }
#
#        # call
#        if verbosity > 0:
#            pass
#        if verbosity > 1:
#            pass
#
#        complete = """"""
#        for key, value in dimension_annotation_strings.items():
#            complete += value
#
#        table_strings[f"{dimension}_table"] += f"""
#    * Number of {dimension}s: {num_dimensions}
#    * {dimension} annotation table features: 
#        {complete}
#        """
#
#    # initialize formatting strings for all enrichment layers
#    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
#    enrichment_strings = {
#        f'{l}':f"""
#    * {l}
#    """
#    for l in enr_layers
#    }
#
#    complete = """"""
#    for key, value in enrichment_strings.items():
#        complete += value
#    table_strings['enrichments'] += complete
#
#    final = """"""
#    for key, value in table_strings.items():
#        final += value
#
#    return final


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
    Query Type                     Condition
    qkey
    sq_1     sample                 Cohort == 2.0
    sq_2     sample  technical_replicate_id > 500
    """

    # TODO make this 'dim', instead of 'Type'
    # TODO run checks, Raise Errors

    # st = ds.sample_table.to_pandas().infer_objects()
    sq = list(query_df.loc[query_df["Type"] == "sample", "Condition"].values)
    sid = sample_id_coordinate_from_query(ds, sq)

    # pt = ds.peptide_table.to_pandas().infer_objects()
    pq = list(query_df.loc[query_df["Type"] == "peptide", "Condition"].values)
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


# TODO dim not table parameter to be consistant
# DEPRECATED
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

    #table = ds[table]
    if table == "sample_table":
        table = get_sample_table(ds)
    else:
        table = get_peptide_table(ds)

    #####?
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
