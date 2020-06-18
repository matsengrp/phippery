"""
@File: PhipData.py

@Author: Jared Galloway

This contains the source code for PhipData
Object and associated functions.
"""

# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import cm

# built-in python3
from functools import reduce
import glob
import os
import re
from collections import defaultdict


def load_counts(
    counts_files,
    peptide_metadata,
    sample_metadata,
    technical_replicate_function="mean",
    add_tech_rep_correlation_to_sample_md=True,
):
    """
    collect data and return PhipData object
    """

    # Load our three tables with helper functions below
    counts, replicate_df = collect_merge_prune_count_data(
        counts_files, technical_replicate_function
    )
    peptide_metadata = collect_peptide_metadata(peptide_metadata)
    sample_metadata = collect_sample_metadata(sample_metadata)

    # this is probably overkill
    assert type(counts) == pd.DataFrame
    assert type(replicate_df) == pd.DataFrame
    assert type(peptide_metadata) == pd.DataFrame
    assert type(sample_metadata) == pd.DataFrame

    # these axis will become xarray coordinates
    assert counts.index.dtype == int
    assert counts.columns.dtype == int
    assert peptide_metadata.index.dtype == int
    assert sample_metadata.index.dtype == int

    # the whole sample metadata could contain extra samples
    # from other references. we only want relevent samples.
    assert set(replicate_df.index).issubset(sample_metadata.index)
    assert set(counts.columns).issubset(sample_metadata.index)
    assert set(counts.columns) == set(replicate_df.index)
    sample_metadata_subset = sample_metadata.loc[list(counts.columns), :]
    assert set(sample_metadata_subset.index) == set(replicate_df.index)
    if add_tech_rep_correlation_to_sample_md:
        sample_metadata = sample_metadata.merge(
            replicate_df, left_index=True, right_index=True
        )

    # TODO we should probably make sure everything lines up
    # TODO this makes the sample and peptide index sorted a requirement.
    # we could simply sort it ourselves here.
    sorted_columns_counts = counts[sorted(counts.columns)]
    assert np.all(sorted_columns_counts.columns == sample_metadata.index)
    assert np.all(sorted_columns_counts.index == peptide_metadata.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    return xr.Dataset(
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
    requirements = ["fastq_pattern", "reference", "experiment"]
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


def collect_merge_prune_count_data(
    counts, technical_replicate_function="mean", pseudo_count_bias=10,
):
    """
    This function takes in a list of paths which
    contains the counts for each peptide alignment
    for each sample. These files should contain
    no header.

    For now, this function only looks for filenames
    which have the pattern
    anything6.3.tsv where 'anything' could be any file
    path or string, the first integer is the sample and
    the second integer is the replicate number.

    :param: counts_dir <str> - the path leading
    to the directory containing counts files
    for each sample and technical replicate

    :param: technical_replicate_function <str>
    - either 'sum' or 'mean'. How we decide to
    summarize both technical replicates for 1 sample.
    """

    technical_replicates = defaultdict(list)
    for f in counts:
        match = re.match(r"\D*(\d+)\.(\d+)\.tsv", os.path.basename(f))
        if match is not None:
            print(
                " ".join(
                    f"processing sample #{match.group(1)}, \
                    technical replicate #{match.group(2)}".split()
                )
            )
            technical_replicates[int(match.group(1))].append(f)

    load = lambda path, sample: pd.read_csv(  # noqa
        path, index_col=0, sep="\t", names=["ID", sample]
    )

    sample_dataframes = {}
    replicate_info = []
    for sample in technical_replicates:
        number_of_technical_replicates = len(technical_replicates[sample])
        if number_of_technical_replicates == 1:
            df = load(technical_replicates[sample][0], sample)
            sample_dataframes[sample] = df
            replicate_info.append([sample, 1, 1.0])
            continue
        elif number_of_technical_replicates == 2:
            rep_1_df = load(technical_replicates[sample][0], sample)
            rep_2_df = load(technical_replicates[sample][1], sample)

            # TODO I don't think they need to have the same indexing order?
            # use set() instead
            assert np.all(set(rep_1_df.index) == set(rep_2_df.index))

            correlation = (
                st.pearsonr(rep_1_df.values.flatten(), rep_2_df.values.flatten())[0]
                if np.any(rep_1_df.values.flatten() != rep_2_df.values.flatten())
                else 1.0
            )
            replicate_info.append([sample, 2, correlation])

            if technical_replicate_function == "sum":
                agg = rep_1_df + rep_2_df
                sample_dataframes[sample] = agg
            elif technical_replicate_function == "mean":
                avg = (rep_1_df + rep_2_df) / 2
                sample_dataframes[sample] = avg
            else:
                raise ValueError(f"{technical_replicate_function} is not implimented")
        else:
            # TODO implement multiple replicates
            # the way to do this is actually to go through all pairs
            # of replicates and compute correlation.
            # you could either take the mean of these or check the threshold
            # for _any_ replicate
            print(
                "warning, greater than two replicates has not been "
                "implemented, skipping sample {sample}, for now"
            )
            continue

    # TODO consider the merge options a little more once we know the
    # nature of the upstream pipeline.
    # are the indexes going to be the same set for every sample?
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        list(sample_dataframes.values()),
    ).fillna(0)

    replicate_df = pd.DataFrame(
        replicate_info, columns=["ID", "num_tech_reps", "tech_rep_correlation"]
    ).set_index("ID")
    assert set(merged_counts_df.columns) == set(replicate_df.index)

    return merged_counts_df, replicate_df
