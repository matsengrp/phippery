"""
@File: PhipData.py

@Author: Jared Galloway

This contains the source code for PhipData
Object and associated functions.
"""

# dependencies
import pandas as pd
import numpy as np
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
    technical_replicate_threshold=0.80,
    technical_replicate_function="mean",
    pseudo_count_bias=10,
    add_tech_rep_correlation_to_sample_md=True,
):
    """
    collect data and return PhipData object
    """

    counts, replicate_df = collect_merge_prune_count_data(
        counts_files,
        technical_replicate_threshold,
        technical_replicate_function,
        pseudo_count_bias,
    )
    peptide_metadata = collect_peptide_metadata(peptide_metadata)
    sample_metadata = collect_sample_metadata(sample_metadata)
    assert set(replicate_df.index).issubset(sample_metadata.index)
    assert set(counts.columns).issubset(sample_metadata.index)
    assert set(counts.columns) == set(replicate_df.index)
    sample_metadata_subset = sample_metadata.loc[list(counts.columns), :]
    assert set(sample_metadata_subset.index) == set(replicate_df.index)
    if add_tech_rep_correlation_to_sample_md:
        sample_metadata = sample_metadata.merge(
            replicate_df, left_index=True, right_index=True
        )

    return PhipData(counts, peptide_metadata, sample_metadata)


def load_csv():
    """
    TODO
    """
    pass


class PhipData(object):
    """
    An object for holding phip-seq datasets.
    """

    def __init__(self, counts, peptide_metadata, sample_metadata):

        self.counts = counts
        self.peptide_metadata = peptide_metadata
        self.sample_metadata = sample_metadata

        self._check_consistancy()

    def _check_consistancy(self):
        """
        check that the inputs have consistant
        indexing, and the required fields.
        """
        assert type(self.counts) == pd.DataFrame
        assert type(self.peptide_metadata) == pd.DataFrame
        assert type(self.sample_metadata) == pd.DataFrame

        assert self.counts.index.dtype == int
        assert self.counts.columns.dtype == int

        # assert that we have metadata for every peptide and sample
        assert (
            len(set(self.counts.columns).difference(set(self.sample_metadata.index)))
            == 0
        )
        assert (
            len(set(self.counts.index).difference(set(self.peptide_metadata.index)))
            == 0
        )

        # TODO where should pruning go?
        # samples that were pruned out
        # pruned_samples = set(sample_metadata.index).difference(set(counts.columns))
        # sample_metadata.drop(pruned_samples, axis=0)

        return None

    def dump_to_csv_stub(self):
        """
        a method which dumps the PhipData
        counts and metadata to a csv file.
        """

        pass

    def subset_stub(self, condition):
        """
        subset the tables based upon some conditional
        string such as

        "experiment == expa"
        """

        pass


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

    # TODO change required field based upon type of experiment. phip or phage-dms
    peptide_metadata = pd.read_csv(peptide_md, sep=",", index_col=0, header=0)
    peptide_metadata.index = peptide_metadata.index.astype(int)
    requirements = ["Oligo"]
    assert np.all([x in peptide_metadata.columns for x in requirements])
    return peptide_metadata


def collect_merge_prune_count_data(
    counts,
    technical_replicate_threshold=0.80,
    technical_replicate_function="mean",
    pseudo_count_bias=10,
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

    :param: technical_replicate_threshold <float>
    - Provide a floating point between -1 and 1
    to use as a correlation threshold before a
    sample is not used

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
