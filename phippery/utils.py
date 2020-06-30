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
import xarray as xr

# built-in python3
import os
import copy


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
    ds, exp_list=[], locus_list=[], locus_name_column="Virus",
):
    """
    subset a phip dataset given a list of experiment and locus names
    """

    if locus_name_column not in ds.peptide_metadata.values:
        raise ValueError(f"{locus_name_column} is not in peptide_metadata coordinate")

    flatten = lambda l: [item for sublist in l for item in sublist]  # noqa

    exp_sample_ids = []
    for exp_name in exp_list:
        exp_sample_ids.append(
            ds.sample_id.where(
                ds.sample_table.loc[:, "experiment"] == exp_name, drop=True
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


def get_all_experiments(ds):
    """ return a list of all available experiments """

    all_exp = ds.sample_table.loc[:, "experiment"]
    return [x for x in set(all_exp.values) if x == x]
