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
    """
    Some tools for PhIP-Seq data analysis. For
    help with a specific command, type phippery COMMAND -h
    """
    pass


# TODO
# xarray option?
# sample/pep meta optional?
@cli.command(
    name="collect-phip-data",
    short_help="Collect and merge counts/metadata into one dataset to be pickle dump'd",
)
@option(
    "-c",
    "--counts_directory",
    required=True,
    type=Path(exists=True),
    help="Directory path to tsv files containing \
        peptide enrichments for each sample. See \
        README.md for directory/file format specifications.\f",
)
@option(
    "-s_meta",
    "--sample_metadata",
    required=True,
    type=Path(exists=True),
    help="File path to sample metadata. See \
    README.md for file format specifications.",
)
@option(
    "-p_meta",
    "--peptide_metadata",
    required=True,
    type=Path(exists=True),
    help="File path to peptide metadata. See \
        README.md for file format specifications.",
)
@option(
    "-o",
    "--output",
    required=True,
    help=" ".join(
        "the file path of the output file \
        where the phip dataset will be pickle \
        dump'd".split()
    ),
)
@option(
    "-tech_rep_thresh",
    "--technical_replicate_correlation_threshold",
    required=False,
    default=0.8,
    show_default=True,
    type=float,
    help=" ".join(
        "This specifies the correlation threshold \
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
    type=str,
    help="This command looks for technical replicates \
        and currently joins them by either summations (option \
        sum) or averaging them (option mean)"
)
@option(
    "-bias",
    "--pseudo_count_bias",
    required=False,
    default=10,
    show_default=True,
    type=int,
    help="This specifies the amount you would like to add"
        "to each peptide count, for each sample, before merging"
        "and aggregating technical replicates"
)
# TODO long description, point to README
def collect_phip_data(
    counts_directory,
    sample_metadata,
    peptide_metadata,
    output,
    technical_replicate_correlation_threshold,
    technical_replicate_function,
    pseudo_count_bias 
):
    """
    """

    # TODO, make a check input format function in utils
    sample_metadata = utils.collect_sample_metadata(sample_metadata)
    peptide_metadata = utils.collect_peptide_metadata(peptide_metadata)

    # TODO obviously these should (maybe not??) be required
    # either that or passed in as an argument.
    # OR, let the user specify a custom name for mockip
    # and library control 
    mock = sample_metadata[sample_metadata["Notes"] == "negative_control"].index
    library_control = sample_metadata[
        sample_metadata["Notes"] == "library_control"
    ].index

    counts = utils.collect_merge_prune_count_data(
        counts_directory,
        technical_replicate_correlation_threshold,
        technical_replicate_function,
        threshold_filter_exceptions=[int(library_control.values), int(mock.values)],
        pseudo_count_bias = pseudo_count_bias
    )

    # xarray Still not sure if this is going to be helpful - I don't think so
    # phip_dataset = utils.create_phip_xarray_dataset(
    #    counts, sample_metadata, peptide_metadata
    # )

    phip_dataset = utils.create_phip_dict_dataset(
        counts, sample_metadata, peptide_metadata
    )

    pickle.dump(phip_dataset, open(output, "wb"))

    return None


# TODO
# * `Group` all analysis into nested command
# for each analysis
#  provide an option for either plotting or dumping
#  analysis data


@cli.command(
    name="fold-analysis",
    short_help="perform standardized fold enrichment analysis on phip dataset",
)
@option(
    "-d",
    "--phip-data",
    required=True,
    type=Path(exists=True),
    help=" ".join(
        "file path to the pickle dump'd phip data object, \
        should have been created by the `collect_phip_data` command \
        ".split()
    ),
)
# TODO, mock/control can be inferred.
@option(
    "-mock",
    "--mock_immunoprecipitation_sample",
    required=True,
    type=int,
    help=" ".join(
        "The beads-only negative control to use when \
        computing the standardized enrichment.".split()
    ),
)
@option(
    "-lib",
    "--library_control_sample",
    required=True,
    type=int,
    help=" ".join(
        "The 'library phage abundance control' to use when \
        computing the standardized enrichment.".split()
    ),
)
@option(
    "-o",
    "--output",
    required=True,
    help=" ".join(
        "the file path of the output file \
        where the standardized enrichment phip \
        dataset will be pickle \
        dump'd".split()
    ),
)
# TODO long description, point to README
def fold_analysis(
    phip_data, mock_immunoprecipitation_sample, library_control_sample, output
):
    """
    """
    import phippery.fold_analysis as fa

    phip_dataset = pickle.load(open(phip_data, "rb"))

    # TODO how to infer the *first* lib control and mock ip.
    # sample_metadata = phip_dataset["sample_metadata"]
    # peptide_metadata = phip_dataset["peptide_metadata"]
    # mock = sample_metadata[
    #    sample_metadata["Notes"] == "negative_control"
    # ].index[0]
    # library_control = sample_metadata[
    #    sample_metadata["Notes"] == "library_control"
    # ].index[0]

    print(f"computing standardized enrichment using {mock_immunoprecipitation_sample}")
    print(f"computing standardized enrichment using {library_control_sample}")

    fa.calculate_fold_normalized_enrichment(
        phip_dataset=phip_dataset,
        library_control_sample=library_control_sample,
        mock_ip_sample=mock_immunoprecipitation_sample,
    )

    pickle.dump(phip_dataset, open(output, "wb"))


# TODO Command plot
# TODO Command analysis
