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
    return None


def stitch_dataset(
    counts, peptide_table, sample_table,
):
    """
    """

    sorted_columns_counts = counts[sorted(counts.columns)]

    # make sure the coordinated match up.
    assert np.all(sorted_columns_counts.columns == sample_table.index)
    assert np.all(sorted_columns_counts.index == peptide_table.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], sorted_columns_counts),
            "sample_table": (["sample_id", "sample_table"], sample_table),
            "peptide_table": (["peptide_id", "peptide_table"], peptide_table),
        },
        coords={
            "sample_id": sorted_columns_counts.columns.values,
            "peptide_id": counts.index.values,
            "sample_table": sample_table.columns,
            "peptide_table": peptide_table.columns,
        },
    )
    pds.attrs["collapsed_variable"] = None
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
    # WE NEED TO MAKE SURE

    load = lambda path, sample: pd.read_csv(  # noqa
        path, index_col=0, sep="\t", names=["sample_id", sample]
    )

    sample_dataframes = [
        load(path, int(os.path.basename(path).split(".")[0])) for path in counts
    ]

    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        sample_dataframes,
    ).fillna(0)

    return merged_counts_df


def add_stats(ds, stats_files):
    """
    add a directory of files describing summary statistics
    for each sample id. tie the infomration into the sample
    table with a column for each row in the summary stats file

    Each summary stat should
    """

    # TODO Finish the docstring here
    # TODO add checks for file format
    # TODO should each contain the same alignment stats? alt, any alignment stat is valid and
    # any sample which doesn't have that get's and NA of sorts  ... ?

    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    alignment_stats = defaultdict(list)
    for sample_alignment_stats in stats_files:
        fp = os.path.basename(sample_alignment_stats)
        sample_id = int(fp.strip().split(".")[0])
        alignment_stats["sample_id"].append(sample_id)
        for line in open(sample_alignment_stats, "r"):
            line = line.strip().split("\t")
            alignment_stats[f"{line[0]}"].append(num(line[1]))

    stats = pd.DataFrame(alignment_stats)
    stats = stats.set_index("sample_id")
    stats = stats.loc[sorted(stats.index)]

    merged = ds.sample_table.combine_first(
        xr.DataArray(stats, dims=["sample_id", "sample_table"])
    )
    return ds.merge(merged)


def dataset_to_csv(ds, file_prefix):
    """
    """

    for dt in list(ds.data_vars):
        ds[f"{dt}"].to_pandas().to_csv(f"{file_prefix}_{dt}.csv", na_rep="NA")


def collect_sample_table(sample_md: str):
    """
    """

    sample_table = pd.read_csv(sample_md, sep=",", index_col=0, header=0)
    sample_table.index = sample_table.index.astype(int)
    requirements = ["fastq_filename", "seq_dir"]
    assert np.all([x in sample_table.columns for x in requirements])
    return sample_table


def collect_peptide_table(peptide_md: str):
    """
    """

    peptide_table = pd.read_csv(peptide_md, sep=",", index_col=0, header=0)
    peptide_table.index = peptide_table.index.astype(int)
    requirements = ["Oligo"]
    assert np.all([x in peptide_table.columns for x in requirements])
    return peptide_table


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
