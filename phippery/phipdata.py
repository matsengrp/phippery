"""
@File: phipdata.py

@Author: Jared Galloway

This contains the source code for PhIP-Seq
Data in the for of an xarray
Object and associated functions for organizing,
collecting, and exporting.
"""


# dependencies
import pandas as pd
import numpy as np
import xarray as xr

# built-in python3
from functools import reduce
import glob
import copy
import os
import re
from collections import defaultdict
import pickle


def load(path):
    """
    simple wrapper for loading xarray datasets from pickle binary
    """

    return pickle.load(open(path, "rb"))


def dump(ds, path):
    """
    simple wrapper for dump'ing xarray datasets to pickle binary
    """

    pickle.dump(ds, open(path, "wb"))


def get_annotation_table(ds, dim="sample"):
    """return a copy of the peptide table after converting all
    the datatypes applying the pandas NaN hueristic"""

    st = (ds[f"{dim}_table"]
        .to_pandas()
        .convert_dtypes()
    )

    return st


def get_sample_table(*args):
    """return a copy of the peptide table after converting all
    the datatypes using the pandas NaN hueristic"""

    return get_annotation_table(*args, dim = "sample")


def get_peptide_table(*args):
    """return a copy of the sample table after converting all
    the datatypes using the pandas NaN heuristic"""

    return get_annotation_table(*args, dim = "peptide")


def dataset_to_wide_csv(ds, file_prefix):
    """
    """
    # TODO make sure we're getting 
    # the right output, compared to input.
    # one ting to do is get rid of the names with data_tables
    # i.e. make another for loop just for the peptide/sample table
    
    
    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    for dt in enr_layers:
        layer = copy.deepcopy(ds[f"{dt}"].to_pandas())
        layer.index.name = ""
        layer.to_csv(
                f"{file_prefix}_{dt}.csv", 
                na_rep="NA", 
                index_label=None
        )

    for at in ["sample", "peptide"]:
        get_annotation_table(ds, dim=at).to_csv(
                f"{file_prefix}_{at}_annotation_table.csv"
        )



def dataset_from_csv(
        counts_matrix_filename, 
        peptide_table_filename, 
        sample_table_filename
):
    """
    """

    counts_df = collect_counts_matrix(counts_matrix_filename)
    peptide_df = collect_peptide_table(peptide_table_filename)
    sample_df = collect_sample_table(sample_table_filename)

    # could do the id coordinate alignment assertion here,
    # as well as the sorting if that's something you'd
    # still like to do (necessary for collapse?)

    return stitch_dataset(
        counts = counts_df, 
        peptide_table=peptide_df, 
        sample_table=sample_df,
    )


def collect_sample_table(sample_table_filename: str):
    """
    """

    sample_table = pd.read_csv(
            sample_table_filename, 
            sep=",",
            index_col=0, 
            header=0
    ) 
    #).convert_dtypes()


    if sample_table.index.name != 'sample_id':
        raise ValueError("The name of the index must be 'sample_id'")

    if sample_table.index.dtype != 'int64':
        raise ValueError("The index values for sample_id must be inferred as integers")


    # TODO let's do away with this requirement
    sample_table.sort_index(inplace=True)

    #requirements = ["fastq_filename", "seq_dir"]
    #assert np.all([x in sample_table.columns for x in requirements])
    return sample_table


def collect_peptide_table(peptide_table_filename: str):
    """
    """
    # TODO 
    # these requirementsare really only for 
    # phip-flow

    # TODO 
    # check on naming for index
    # i.e. peptide_table.index.name = "peptide_id"

    peptide_table = pd.read_csv(
            peptide_table_filename, 
            sep=",", 
            index_col=0, 
            header=0
    )

    if peptide_table.index.name != 'peptide_id':
        raise ValueError

    if peptide_table.index.dtype != 'int64':
        raise ValueError("The index values for peptide_id must be inferred as integers")

    peptide_table.sort_index(inplace=True)
    return peptide_table


def collect_counts_matrix(counts_matrix_filename: str):
    """
    """
    # TODO should we throw a warning if not integer?

    counts_matrix = pd.read_csv(
        counts_matrix_filename, 
        sep=",", 
        index_col=0, 
        header=0
    )

    # TODO check matrix type is all numeric
    # and/or do another try except block to force them?

    try:
        counts_matrix.columns = counts_matrix.columns.astype(int)
    except:
        raise ValueError("column header values much be able to cast to type 'int' to match peptide table index")

    try:
        counts_matrix.index = counts_matrix.index.astype(int)
    except:
        raise ValueError("row index values much be able to cast to type 'int' to match peptide table index")

    counts_matrix.sort_index(inplace=True)
    counts_matrix.sort_index(axis=1, inplace=True)

    return counts_matrix
