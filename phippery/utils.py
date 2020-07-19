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


def get_std_dev_non_control(
    ds, control_status_column="control_status", sample_label="empirical"
):
    """
    given a dataset, compute the standard deviation on all biological, non-control samples
    defined by the 'control_status_column' and the correct factor, 'sample_label'.
    """

    if control_status_column not in ds.sample_metadata.values:
        raise ValueError(
            f"{control_status_column} is not in sample_metadata coordinate"
        )

    # cut controls
    empirical_sample_ids = ds.sample_id.where(
        ds.sample_table.loc[:, control_status_column] == sample_label, drop=True
    ).values

    empirical_ds = ds.loc[dict(sample_id=empirical_sample_ids)]

    return np.std(empirical_ds.counts.values)


def get_bucket_index(buckets, v):

    for i in range(len(buckets) - 1):
        if v in range(buckets[i], buckets[i + 1]):
            return i

    return None


def get_exp_prot_subset(
    ds,
    exp_list=[],
    locus_list=[],
    exp_name_column="experiment",
    locus_name_column="Virus",
):
    """
    subset a phip dataset given a list of experiment and locus names
    """
    # TODO
    print("Deprecated")

    if locus_name_column not in ds.peptide_metadata.values:
        raise ValueError(f"{locus_name_column} is not in peptide_metadata coordinate")

    flatten = lambda l: [item for sublist in l for item in sublist]  # noqa

    exp_sample_ids = []
    for exp_name in exp_list:
        exp_sample_ids.append(
            ds.sample_id.where(
                ds.sample_table.loc[:, exp_name_column] == exp_name, drop=True
            ).values
        )
    exp_sample_ids = flatten(exp_sample_ids)
    if len(exp_sample_ids) == 0:
        exp_sample_ids = list(ds.sample_id.values)

    locus_peptide_ids = []
    for locus_name in locus_list:
        locus_peptide_ids.append(
            ds.peptide_id.where(
                ds.peptide_table.loc[:, locus_name_column] == locus_name, drop=True
            ).values
        )
    locus_peptide_ids = flatten(locus_peptide_ids)
    if len(locus_peptide_ids) == 0:
        locus_peptide_ids = list(ds.peptide_id.values)

    ds_copy = copy.deepcopy(ds)

    return ds_copy.loc[dict(sample_id=exp_sample_ids, peptide_id=locus_peptide_ids)]


def get_all_sample_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a sample table column
    """

    all_exp = ds.sample_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def get_all_peptide_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a peptide table column
    """

    all_exp = ds.peptide_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def pairwise_correlation_by_sample_group(ds, group="sample_ID"):
    """
    a method which computes pairwise cc for all
    sample in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of samples
    in the group.
    """

    if group not in ds.sample_metadata.values:
        raise ValueError("{group} does not exist in sample metadata")

    groups, pw_cc, n = [], [], []
    for group, group_ds in ds.groupby(ds.sample_table.loc[:, group]):
        groups.append(group)
        n.append(len(group_ds.sample_id.values))
        if len(group_ds.sample_id.values) < 2:
            pw_cc.append(0)
            continue
        correlations = []
        for sample_ids in itertools.combinations(group_ds.sample_id.values, 2):
            sample_0_enrichment = group_ds.counts.loc[:, sample_ids[0]]
            sample_1_enrichment = group_ds.counts.loc[:, sample_ids[1]]
            correlation = (
                st.pearsonr(sample_0_enrichment, sample_1_enrichment)[0]
                if np.any(sample_0_enrichment != sample_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(round(sum(correlations) / len(correlations), 3))

    return pd.DataFrame({"group": groups, "pairwise_correlation": pw_cc, "n_reps": n})


def melted_mean_collapse_groups(ds, by="sample_ID"):
    """
    This function melts a dataset after taking the average of
    all values in sample metadata coordinate group.
    The nature of this function forces us to throw out sample
    table information for all but the first of the group
    encountered in the dataset - defined by the "by" category.
    """

    # create a "by" sample table, which only takes the first instance
    # of a group in the sample metadata, and sets that as the index.
    sample_table = ds.sample_table.to_pandas()
    to_drop = []
    bio_id = []
    for sam_id, row in sample_table.iterrows():
        if row[by] in bio_id:
            to_drop.append(sam_id)
        bio_id.append(row[by])
    bio_sample_table = sample_table.drop(to_drop).set_index(by, verify_integrity=True)

    # a new counts table, melty af.
    avg_bio_reps_df = ds.counts.groupby(ds.sample_table.loc[:, by]).mean().to_pandas()
    melted = avg_bio_reps_df.reset_index().melt(id_vars=["peptide_id"])
    melted.rename({"sample_table": by}, axis=1, inplace=True)

    # just our peptide table - it's perfect how it is except we need the pid to merge by.
    peptide_table = ds.peptide_table.to_pandas().reset_index()

    melted_peptide = melted.merge(peptide_table, on="peptide_id")
    all_melted = melted_peptide.merge(bio_sample_table, on=by)

    return all_melted


def melt_ds(ds):
    """
    return an new dataset, compltely melted. Ideal for ggplotting.
    """
    counts = ds.counts.to_pandas().reset_index().melt(id_vars=["peptide_id"])
    peptide_table = ds.peptide_table.to_pandas().reset_index()
    sample_table = ds.sample_table.to_pandas().reset_index()
    counts_peptide = counts.merge(peptide_table, on="peptide_id")
    return counts_peptide.merge(sample_table, on="sample_id")


def iter_sample_groups(ds, column):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the sample metadata coodinate.
    """
    for group, group_ds_idx in ds.sample_id.groupby(ds.sample_table.loc[:, column]):
        group_ds = ds.loc[dict(sample_id=group_ds_idx.sample_id.values)]
        yield group, group_ds


def iter_peptide_groups(ds, column):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the peptide metadata coodinate.
    """
    for group, group_ds_idx in ds.peptide_id.groupby(ds.peptide_table.loc[:, column]):
        group_ds = ds.loc[dict(peptide_id=group_ds_idx.peptide_id.values)]
        yield group, group_ds


def sample_id_subset(ds, where, is_equal_to):
    """
    grab the subset of sample id's representing a factor
    group in the sample metadata
    """
    return list(
        ds.sample_id.where(
            ds.sample_table.loc[:, where] == is_equal_to, drop=True
        ).sample_id.values
    )


def peptide_id_subset(ds, where, is_equal_to):
    """
    grab the subset of sample id's representing a factor
    group in the peptide metadata
    """
    return list(
        ds.peptide_id.where(
            ds.peptide_table.loc[:, where] == is_equal_to, drop=True
        ).peptide_id.values
    )


def subset_ds_by_sample_factor(ds, where, is_equal_to):
    """
    grab the subset loc of the ds representing a factor
    group in the sample metadata
    """
    sample_ids = list(
        ds.sample_id.where(
            ds.sample_table.loc[:, where] == is_equal_to, drop=True
        ).sample_id.values
    )

    return ds.loc[dict(sample_id=sample_ids)]


def subset_ds_by_peptide_factor(ds, where, is_equal_to):
    """
    grab the subset loc of the ds representing a factor
    group in the peptide metadata
    """
    peptide_ids = list(
        ds.peptide_id.where(
            ds.peptide_table.loc[:, where] == is_equal_to, drop=True
        ).peptide_id.values
    )

    return ds.loc[dict(peptide_id=peptide_ids)]
