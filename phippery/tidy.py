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


def tidy_ds(ds, sample_metadata_columns="all", peptide_metadata_columns="all"):
    """
    return an new dataset, with all data compltely melted.
    This means that instead of a counts matrix
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
    if sample_metadata_columns != "all":
        sample_table = sample_table[["sample_id"] + sample_metadata_columns]
    if peptide_metadata_columns != "all":
        peptide_table = peptide_table[["peptide_id"] + peptide_metadata_columns]

    # merge the metadata into the melted datatables
    data_peptide = merged_counts_df.merge(peptide_table, on="peptide_id")
    return data_peptide.merge(sample_table, on="sample_id")
