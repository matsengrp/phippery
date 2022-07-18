"""
@File: utils.py

@Author: Jared Galloway
"""

# local
from phippery.utils import iter_sample_groups

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


def tidy_ds(ds):
    """
    return the provided dataset in 'tall form'.

    This means that all information from the entire dataset and
    all the tables it includes will be 'melted' into a single dataframe.
    This format is far less effecient in terms of storage, and should not
    be used on the entire dataset, but rather, some subset of the dataset
    when you are ready to plot the relevent information.

    The means that for each row, you will get the following columns:

    sample_id,
    peptide_id,
    all sample metadata (a column for each),
    all peptide metadata (a column for each),
    all enrichment tables (a column for each).

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
    # sample_coord_dim = ds.attrs["sample_coord_dim"]
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, on=["peptide_id", "sample_id"]), melted_data
    )

    # grab only the columns of the metadata we want in the resulting dataframe
    peptide_table = ds.peptide_table.to_pandas().reset_index().infer_objects()
    sample_table = ds.sample_table.to_pandas().reset_index().infer_objects()

    # merge the metadata into the melted datatables
    data_peptide = merged_counts_df.merge(peptide_table, on="peptide_id")
    return data_peptide.merge(sample_table, on="sample_id")
