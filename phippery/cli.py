"""
@File: phippery_cli.py

@Author: Jared Galloway

Command line interface (CLI) for phippery.
"""

# built-in
import pickle

# dependencies
import pandas as pd
import numpy as np
from click import Choice, Path, command, group, option
import click

# local
import phippery.utils as utils
import phippery.phipdata as phipdata

# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Some tools for PhIP-Seq data analysis. For
    help with a specific command, type phippery COMMAND -h
    """
    pass


@cli.command(name="collect-phip-data")
@click.argument("counts", required=True, nargs=-1, type=Path(exists=True))
@option(
    "-s_meta",
    "--sample_metadata",
    required=True,
    type=Path(exists=True),
    help="File path to sample metadata. See README.md for file format specifications.",
)
@option(
    "-p_meta",
    "--peptide_metadata",
    required=True,
    type=Path(exists=True),
    help="File path to peptide metadata. See README.md for file format specifications.",
)
@option(
    "-o",
    "--output",
    required=True,
    help="the file path of the output file where the phip dataset will be pickle dump'd",
)
def collect_phip_data(
    counts, sample_metadata, peptide_metadata, output,
):
    """
    Collect sample and peptide metadata tables along with individual two-column tsv files, for each sample, and produce a properly formatted xarray dataset.
    """

    print(list(counts))
    xds = phipdata.counts_metadata_to_dataset(
        counts_files=list(counts),
        peptide_metadata=open(peptide_metadata, "r"),
        sample_metadata=open(sample_metadata, "r"),
    )

    pickle.dump(xds, open(output, "wb"))

    return None


@cli.command()
@option(
    "-d",
    "--peptide-metadata",
    type=Path(exists=True),
    required=True,
    help="The file path to the peptide metadata",
)
@option(
    "-o",
    "--output-fasta",
    type=Path(),
    required=True,
    help="The output file path for the fasta",
)
def peptide_md_to_fasta(peptide_metadata, output_fasta):
    """
    convert peptide metadata to fasta format.

    For each peptide, we will add an entry to a fasta file with
    the unique peptide_id as the only entry into '>' header and
    oligonucleotide encoding on the line below.
    """
    utils.convert_peptide_metadata_to_fasta(peptide_metadata, output_fasta)
