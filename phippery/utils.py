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
import numpy as np
import scipy.stats as st

# built-in python3
from functools import reduce
import glob
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import cm


def create_phip_dict_dataset(
    counts: pd.DataFrame, sample_metadata: pd.DataFrame, peptide_metadata: pd.DataFrame
):
    """
    simply drop samples not included in the counts,
    assert the indexing is consistant,
    """
    print(counts.index)

    assert counts.index.dtype == int
    assert counts.columns.dtype == int

    # assert that we have metadata for every peptide and sample
    assert len(set(counts.columns).difference(set(sample_metadata.index))) == 0
    assert len(set(counts.index).difference(set(peptide_metadata.index))) == 0

    # samples that were pruned out
    pruned_samples = set(sample_metadata.index).difference(set(counts.columns))
    sample_metadata.drop(pruned_samples, axis=0)

    # should we do this for peptides too?

    return {
        "counts": counts,
        "sample_metadata": sample_metadata,
        "peptide_metadata": peptide_metadata,
    }


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


# TODO Is there a scenerio where this functionality
# would be nice?
def collect_merge_count_data():
    """ dont throw out technical reps? """
    pass


def collect_merge_prune_count_data(
    counts,
    technical_replicate_threshold=0.80,
    technical_replicate_function="mean",
    threshold_filter_exceptions=[],
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
    summarize both technical replicates for 1 sample."""

    technical_replicates = defaultdict(list)
    # for f in glob.iglob(os.path.join(counts_dir, "*.tsv")):
    print(counts)
    for f in counts:
        print(f)
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

            # TODO I don't think they need to have the same indexing order?
            assert np.all(rep_1_df.index == rep_2_df.index)

            # TODO should we store all the samples correlation for plotting?
            # or at least make it an option
            if (
                st.pearsonr(rep_1_df.values.flatten(), rep_2_df.values.flatten())[0]
                < technical_replicate_threshold
                and sample not in threshold_filter_exceptions
            ):
                print(
                    f" ".join(
                        f"throwing out sample: {sample} for having \
                        a technical replicate correlation lower than \
                        {technical_replicate_threshold}".split()
                    )
                )
                continue

            rep_1_df += pseudo_count_bias
            rep_2_df += pseudo_count_bias
            agg = rep_1_df + rep_2_df
            sample_dataframes[sample] = (
                agg if technical_replicate_function == "sum" else agg / 2
            )
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

    return merged_counts_df


def plot_peptide_enrichment_by_nt_position(
    ds, strain_pattern, sample, out, cmap="inferno"
):
    """
    This function takes plots a subset of tile enrichments
    as a function of nucleotide position (as described by
    nt_start, nt_end in peptide metadata) for all specific sample.

    # TODO bad practice to have a variable take on multiple types.
    :param: strain_pattern -
    if this parameter is a string, it is assumed to be a regex pattern
    and all virus strains in the ds that match are plotted. Otherwise,
    a list of viruses should be provided to
    """
    # TODO assert phip dataset consistency

    # if type(strain_pattern) == str:
    all_strains = set(ds["peptide_metadata"]["Virus_Strain"])
    selected_strains = [
        re.match(strain_pattern, f"{strain}")[0]
        for strain in all_strains
        if re.match(strain_pattern, f"{strain}") is not None
    ]
    assert len(all_strains) != 0

    fig, ax = plt.subplots()
    color_map = cm.get_cmap(cmap, len(selected_strains)).colors
    for i, strain in enumerate(selected_strains):

        strain_indices = ds["peptide_metadata"][
            ds["peptide_metadata"]["Virus_Strain"] == strain
        ].index

        tile_start = ds["peptide_metadata"]["nt_start"][strain_indices]
        tile_end = ds["peptide_metadata"]["nt_end"][strain_indices]
        enrichment = ds["counts"][sample][strain_indices]
        lines = [
            [(x1, y), (x2, y)] for x1, x2, y in zip(tile_start, tile_end, enrichment)
        ]
        lc = mc.LineCollection(lines, linewidths=1, label=strain, color=color_map[i])
        ax.add_collection(lc)

    ax.autoscale()
    sample_info = ds["sample_metadata"]["sample_info"][sample]
    seroepi_paper_id = ds["sample_metadata"]["seroepi_paper_id"][sample]
    notes = ds["sample_metadata"]["Notes"][sample]
    ax.set_title(
        f"sample_id : {sample},\n \
        sample_info : {sample_info},\n \
        seroepi_paper_id : {seroepi_paper_id},\n \
        Notes : {notes}"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_ylabel(f"Standardized fold enrichment")
    ax.set_xlabel(f"peptide tile, in order (genome position)")
    # ax.grid()
    plt.tight_layout()
    fig.savefig(f"{out}")
    plt.close(fig)


# TODO implement
def check_phip_dataset_consistancy_stub(phip_dataset):
    """
    do some error checking on an xarray
    phip dataset. This should probably consist of

    1. Make sure the dimensions are correct
    2. Check that we have metadata for every sample and peptide
    3. ???? TODO

    """
    data_consistancy = True
    return data_consistancy


def convert_peptide_metadata_to_fasta(peptide_metadata, out):
    """
    Take in peptide metadata dataframe, and write a fasta
    format representation of the oligos
    """

    def trim_index(sequence):
        return "".join([nt for nt in sequence if nt.isupper()])

    fasta_fp = open(out, "w")
    peptide_metadata = pd.read_csv(peptide_metadata, index_col=0, header=0)
    requirements = ["Oligo"]
    assert peptide_metadata.index.name == "ID"
    assert np.all([x in peptide_metadata.columns for x in requirements])
    for index, row in peptide_metadata.iterrows():
        ref_sequence = trim_index(row["Oligo"])
        fasta_fp.write(f">{index}\n{ref_sequence}\n")
