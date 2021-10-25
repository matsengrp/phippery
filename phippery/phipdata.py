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
        .infer_objects()
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


def stitch_dataset(
    counts, peptide_table, sample_table,
):
    """given all corrected datatypes"""

    # TODO Today - this may cause issues, why do we need to sort, again?
    # and if we do need to sort, can we simply sort the matrix, as well?
    # sorted_columns_counts = counts[sorted(counts.columns)]

    # make sure the coordinated match up.
    assert np.all(counts.columns == sample_table.index)
    assert np.all(counts.index == peptide_table.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], counts),
            "sample_table": (["sample_id", "sample_metadata"], sample_table),
            "peptide_table": (["peptide_id", "peptide_metadata"], peptide_table),
        },
        coords={
            "sample_id": counts.columns.values,
            "peptide_id": counts.index.values,
            "sample_metadata": sample_table.columns.values,
            "peptide_metadata": peptide_table.columns.values,
        },
    )
    return pds


def collect_merge_prune_count_data(counts):
    """
    This function takes in a list of paths which
    contains the counts for each peptide alignment
    for each sample. These files should contain
    no header.

    :param: counts <str> - a list of paths leading
    to raw peptide enrichment counts for each sample
    """

    # TODO ADD CHECKS
    # WE NEED TO MAKE SURE EACH FOLLOWS A CERTAIN FORMAT
    
    # TODO remove prune from the name.

    load = lambda path, sample: pd.read_csv(  # noqa
        path, index_col=0, sep="\t", names=["sample_id", sample]
    )

    sample_dataframes = [
        load(path, int(os.path.basename(path).split(".")[0])) for path in counts
    ]

    # TODO do we really need to fill na with 0?
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        sample_dataframes,
    ).fillna(0)

    merged_counts_df.columns = merged_counts_df.columns.astype(int)
    merged_counts_df.index = merged_counts_df.index.astype(int)
    merged_counts_df.sort_index(inplace=True)
    merged_counts_df.sort_index(axis=1, inplace=True)

    return merged_counts_df


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
        print(layer.index.name)
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
    ).convert_dtypes()
    #) 


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
    ).convert_dtypes()
    #).convert_dtypes()

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
    ).convert_dtypes()

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


def convert_peptide_table_to_fasta(peptide_table, out):
    """
    Take in peptide metadata dataframe, and write a fasta
    format representation of the oligos
    """

    fasta_fp = open(out, "w")
    peptide_table = pd.read_csv(peptide_table, index_col=0, header=0)
    requirements = ["Oligo"]
    assert peptide_table.index.name == "peptide_id"
    assert np.all([x in peptide_table.columns for x in requirements])
    for index, row in peptide_table.iterrows():
        ref_sequence = trim_index(row["Oligo"])
        fasta_fp.write(f">{index}\n{ref_sequence}\n")


def trim_index(sequence):
    return "".join([nt for nt in sequence if nt.isupper()])
