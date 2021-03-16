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
import scipy.stats as st

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


def add_stats(ds, stats_files):
    """
    add a directory of files describing summary statistics
    for each sample id. tie the infomration into the sample
    table with a column for each row in the summary stats file

    Each summary stat should
    """

    # TODO add checks

    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    alignment_stats = defaultdict(list)
    # for sample_alignment_stats in glob.glob(file_pattern):
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
        xr.DataArray(stats, dims=["sample_id", "sample_metadata"])
    )
    return ds.merge(merged)


def counts_metadata_to_dataset(
    counts_files, peptide_metadata, sample_metadata,
):
    """
    collect data and return an organized xarray object.
    This function assumes the counts are all seperate tsv files
    and uses the function 'collect_merge_prune_count_data' below
    to collect them into a dataframe, finally merging them all into
    a single dataset.
    """

    counts = collect_merge_prune_count_data(counts_files)
    peptide_metadata = collect_peptide_metadata(peptide_metadata)
    sample_metadata = collect_sample_metadata(sample_metadata)

    assert set(counts.columns).issubset(sample_metadata.index)
    sample_metadata = sample_metadata.loc[sorted(counts.columns), :]
    sorted_columns_counts = counts[sorted(counts.columns)]
    assert np.all(sorted_columns_counts.columns == sample_metadata.index.values)
    assert np.all(sorted_columns_counts.index == peptide_metadata.index.values)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], sorted_columns_counts),
            "sample_table": (["sample_id", "sample_metadata"], sample_metadata),
            "peptide_table": (["peptide_id", "peptide_metadata"], peptide_metadata),
        },
        coords={
            "sample_id": sorted_columns_counts.columns.values,
            "peptide_id": counts.index.values,
            "sample_metadata": sample_metadata.columns,
            "peptide_metadata": peptide_metadata.columns,
        },
    )
    pds.attrs["collapsed_variable"] = None
    return pds


def dataset_to_csv(ds, file_prefix):
    """
    a method to dump the relevent tables to csv,
    given an `xarray.dataset` containing phip data.
    """
    for dt in list(ds.data_vars):
        ds[f"{dt}"].to_pandas().to_csv(f"{file_prefix}_{dt}.csv", na_rep="NA")


def csv_to_dataset(
    counts, peptide_table, sample_table,
):
    """
    collect data and return and `xarray.dataset` object,
    organized by respective corrdinates.

    This is similar to the function 'counts_metadata_to_dataset'
    but takes in csv's for all of the tables, assuming the
    counts table has already been generated, and the indices match
    the repective corrdinated for the metadata tables.
    """

    # Load our three tables with helper functions below
    counts = pd.read_csv(counts, sep=",", index_col=0, header=0)
    counts.index = counts.index.astype(int)
    counts.columns = counts.columns.astype(int)
    peptide_metadata = collect_peptide_metadata(peptide_table)
    sample_metadata = collect_sample_metadata(sample_table)

    # these axis will become xarray coordinates
    assert peptide_metadata.index.dtype == int
    assert sample_metadata.index.dtype == int

    sorted_columns_counts = counts[sorted(counts.columns)]
    assert np.all(sorted_columns_counts.columns == sample_metadata.index)
    assert np.all(sorted_columns_counts.index == peptide_metadata.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], sorted_columns_counts),
            "sample_table": (["sample_id", "sample_metadata"], sample_metadata),
            "peptide_table": (["peptide_id", "peptide_metadata"], peptide_metadata),
        },
        coords={
            "sample_id": sorted_columns_counts.columns.values,
            "peptide_id": counts.index.values,
            "sample_metadata": sample_metadata.columns,
            "peptide_metadata": peptide_metadata.columns,
        },
    )
    pds.attrs["collapsed_variable"] = None
    return pds


def df_to_dataset(
    counts_df, peptide_table_df, sample_table_df,
):

    # these axis will become xarray coordinates
    assert counts_df.index.dtype == int
    assert counts_df.columns.dtype == int
    assert peptide_table_df.index.dtype == int
    assert sample_table_df.index.dtype == int

    sorted_columns_counts_df = counts_df[sorted(counts_df.columns)]
    assert np.all(sorted_columns_counts_df.columns == sample_table_df.index)
    assert np.all(sorted_columns_counts_df.index == peptide_table_df.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], sorted_columns_counts_df),
            "sample_table": (["sample_id", "sample_metadata"], sample_table_df),
            "peptide_table": (["peptide_id", "peptide_metadata"], peptide_table_df),
        },
        coords={
            "sample_id": sorted_columns_counts_df.columns.values,
            "peptide_id": counts_df.index.values,
            "sample_metadata": sample_table_df.columns,
            "peptide_metadata": peptide_table_df.columns,
        },
    )
    pds.attrs["collapsed_variable"] = None
    return pds


def collect_sample_metadata(sample_md: str):
    """
    simply load in the sample metadata
    and return the xarray DataArray with correct
    dimensions.

    This could certainly be extended in the future.
    Mainy for checking data format consistancy?
    """

    sample_metadata = pd.read_csv(sample_md, sep=",", index_col=0, header=0)
    sample_metadata.index = sample_metadata.index.astype(int)
    requirements = ["fastq_filename", "reference", "seq_dir"]
    assert np.all([x in sample_metadata.columns for x in requirements])
    return sample_metadata


def collect_peptide_metadata(peptide_md: str):
    """
    simply load in the peptide metadata
    and return the pandas array with correct
    dimensions.

    This could certainly be extended in the future.
    """

    peptide_metadata = pd.read_csv(peptide_md, sep=",", index_col=0, header=0)
    peptide_metadata.index = peptide_metadata.index.astype(int)
    requirements = ["Oligo"]
    assert np.all([x in peptide_metadata.columns for x in requirements])
    return peptide_metadata


def collect_merge_prune_count_data(counts):
    """
    This function takes in a list of paths which
    contains the counts for each peptide alignment
    for each sample. These files should contain
    no header.

    :param: counts <str> - a list of paths leading
    to raw peptide enrichment counts for each sample
    """

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
