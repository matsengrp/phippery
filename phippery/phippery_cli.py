"""
@File: phippery_cli.py

@Author: Jared Galloway

Command line interface (CLI) for phippery.
"""

# built-in
import pickle

# dependencies
import pandas as pd
import xarray as xr
import numpy as np
from click import Choice, Path, command, group, option

# local
import phippery.utils as utils

# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """some tools for phip-seq analysis"""
    pass


@cli.command(name="collect_phip_data")
@option(
    "-c",
    "--counts_directory",
    required=True,
    help="Directory path to tsv files containing \
        peptide enrichments for each ",
)
@option(
    "-s_meta", "--sample_metadata", required=True, help="File path to sample metadata",
)
@option(
    "-p_meta",
    "--peptide_metadata",
    required=True,
    help="File path to peptide metadata",
)
@option(
    "-o",
    "--output",
    required=True,
    help="the name of the output file to dump the dataset",
)
@option(
    "-tech_rep_thresh",
    "--technical_replicate_correlation_threshold",
    required=False,
    default=0.8,
    show_default=True,
    type=float,
    help=" ".join(
        "This species the the correlation threshold \
            which must be met before simply not including the \
            sample in the analysis dataset.".split(),
    ),
)
@option(
    "-tech_rep_agg",
    "--technical_replicate_function",
    required=False,
    default="mean",
    show_default=True,
    help=" ".join(
        "This command looks for technical replicates \
            and currently joins them by euther summations (option \
            'sum') or averaging them (option mean)".split(),
    ),
)
def collect_phip_data(
    counts_directory,
    sample_metadata,
    peptide_metadata,
    output,
    technical_replicate_correlation_threshold,
    technical_replicate_function,
):
    """collect the dataset and dump it to output"""

    # TODO, make a check input format function in utils
    counts_xr = utils.collect_merge_prune_count_data(
        counts_directory,
        technical_replicate_correlation_threshold,
        technical_replicate_function,
    )
    sample_metadata_xr = utils.collect_sample_metadata(sample_metadata)
    peptide_metadata_xr = utils.collect_peptide_metadata(peptide_metadata)
    phip_dataset = utils.create_phip_dataset(
        counts_xr, sample_metadata_xr, peptide_metadata_xr
    )
    pickle.dump(phip_dataset, open(output, "wb"))

    return None
