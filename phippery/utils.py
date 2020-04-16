"""
@File: utils.py

@Author: Jared Galloway

UNDER CONSTRUCTION

This file will include some helpful functions
for the phippery package CLI. The primary
data struct we put counts data into is

So far it includes functions to"

* compile tsv files into a phip dataset
* TODO check phip_dataset attributed

"""

# dependencies
import pandas as pd
import xarray as xr
import numpy as np
import scipy.stats as st

# built-in python3
from functools import reduce
import glob
import os
import re
from collections import defaultdict


def create_phip_dataset(
    counts: xr.DataArray, sample_metadata: xr.DataArray, peptide_metadata: xr.DataArray
):
    """
    Merge the three tsv into one xr DataSet
    after doing some checks on the indexing.
    """

    # TODO always make sure we have the same
    # data types going on
    assert np.all(counts.peptide_id == peptide_metadata.peptide_id)
    assert np.all(counts.sample_id == sample_metadata.sample_id)

    # TODO add some attributes
    phip_dataset = xr.Dataset(
        {"counts": counts, "peptides": peptide_metadata, "samples": sample_metadata}
    )

    return phip_dataset


def collect_sample_metadata(sample_md: str):
    """
    simply load in the sample metadata
    and return the xarray DataArray with correct
    dimensions.

    This could certainly be extended in the future.
    """
    sample_metadata = pd.read_csv(sample_md, sep="\t", index_col=0, header=0)
    return xr.DataArray(sample_metadata, dims=["sample_id", "sample_metadata"])


def collect_peptide_metadata(peptide_md: str):
    """
    simply load in the peptide metadata
    and return the xarray DataArray with correct
    dimensions.

    This could certainly be extended in the future.
    """
    peptide_metadata = pd.read_csv(peptide_md, sep="\t", index_col=0, header=0)
    return xr.DataArray(peptide_metadata, dims=["peptide_id", "peptide_metadata"])


def collect_merge_prune_count_data(
    counts_dir: str,
    technical_replicate_threshold=0.80,
    technical_replicate_function="mean",
):
    """
    This function takes in a directory path which
    contains the counts for each peptide alignment
    for each sample. These files should contain
    no header and

    For now, this function only looks for filenames
    which have the pattern
    r'\D*(\d+)\.(\d+)\.tsv' # noqa
    where the first integer is the sample and
    the second integer is the replicate number.

    :param: counts_dir <str> - the path leading
    to the directory containing counts files
    for each sample and technical replicate

    :param: technical_replicate_threshold <float>
    - Provide a floating point between -1 and 1
    to use as a correlation threshold before a
    sample is not used

    :param: tech
    """

    technical_replicates = defaultdict(list)
    for f in glob.iglob(os.path.join(counts_dir, "*.tsv")):
        match = re.fullmatch(r"\D*(\d+)\.(\d+)\.tsv", f)
        if match is not None:
            print(
                " ".join(
                    f"processing sample #{match.group(1)}, \
                    technical replicate #{match.group(2)}".split()
                )
            )
            technical_replicates[int(match.group(1))].append(f)

    load = lambda path, sample: pd.read_csv(  # noqa
        path, index_col=0, sep="\t", names=["id", sample]
    )

    sample_dataframes = {}
    for sample in technical_replicates:
        number_of_technical_replicates = len(technical_replicates[sample])
        if number_of_technical_replicates == 1:
            sample_dataframes[sample] = load(technical_replicates[sample][0], sample)
            continue
        elif number_of_technical_replicates == 2:
            rep_1_df = load(technical_replicates[sample][0], sample)
            rep_2_df = load(technical_replicates[sample][1], sample)

            # TODO I dont think they need to have the same indexing order?
            assert np.all(rep_1_df.index == rep_2_df.index)

            if (
                st.pearsonr(rep_1_df.values.flatten(), rep_2_df.values.flatten())[0]
                < technical_replicate_threshold
            ):
                continue

            agg = rep_1_df + rep_2_df
            sample_dataframes[sample] = (
                agg if technical_replicate_function == "sum" else agg / 2
            )
        else:
            # TODO implement multiple replicates
            print(
                " ".join(
                    f"warning, greater than two replicates \
                    has not been implemented, skipping sample \
                    {sample}, for now".split()
                )
            )
            continue

    # TODO consider the merge options a little more once we know the
    # nature of the upstream pipeline.
    # are the indexes going to be the same set for every sample?
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        list(sample_dataframes.values()),
    ).fillna(0)

    return xr.DataArray(merged_counts_df, dims=["peptide_id", "sample_id"])


# TODO implement
def check_phip_dataset_consistancy_stub(phip_dataset: xr.Dataset):
    """
    do some error checking on an xarray
    phip dataset. This should probably consist of

    1. Make sure the dimensions are correct
    2. Check that we have metadata for every sample and peptide
    3. ???? TODO

    """
    data_consistancy = True
    return data_consistancy
