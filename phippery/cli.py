"""
@File: phippery_cli.py

@Author: Jared Galloway

Command line interface (CLI) for phippery.
"""

# built-in
import pickle
import glob

# dependencies
import pandas as pd
import numpy as np
from click import Choice, Path, command, group, option
import click

# local
from phippery.phipdata import convert_peptide_table_to_fasta
from phippery.phipdata import collect_merge_prune_count_data
from phippery.phipdata import add_stats
from phippery.phipdata import collect_sample_table
from phippery.phipdata import collect_peptide_table
from phippery.phipdata import stitch_dataset
from phippery.phipdata import dump


# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Some tools for PhIP-Seq data analysis. For
    help with a specific command, type phippery COMMAND -h
    """
    pass


@cli.command(name="collect-phip-data")
@option(
    "-s_meta",
    "--sample_table",
    required=True,
    type=Path(exists=True),
    help="File path to sample metadata. See README.md for file format specifications.",
)
@option(
    "-p_meta",
    "--peptide_table",
    required=True,
    type=Path(exists=True),
    help="File path to peptide metadata. See README.md for file format specifications.",
)
@option(
    "-c",
    "--counts-file-pattern",
    required=True,
    help="File pattern (glob) for all counts files to be merged",
)
@option(
    "-s",
    "--stats-file-pattern",
    required=True,
    help="File pattern (glob) for all counts files to be merged",
)
@option(
    "-o",
    "--output",
    required=True,
    help="the file path of the output file where the phip dataset will be pickle dump'd",
)
def collect_phip_data(
    sample_table, peptide_table, counts_file_pattern, stats_file_pattern, output,
):
    """
    Collect sample and peptide metadata tables along with individual two-column tsv files, for each sample, and produce a properly formatted xarray dataset.
    """

    counts = [f for f in glob.glob(counts_file_pattern)]
    stats = [f for f in glob.glob(stats_file_pattern)]

    merged_counts = collect_merge_prune_count_data(counts)
    peptide_table = collect_peptide_table(peptide_table)
    sample_table = collect_sample_table(sample_table)

    ds = stitch_dataset(
        counts_files=merged_counts,
        peptide_table=peptide_table,
        sample_table=sample_table,
    )

    ds = add_stats(ds, stats)
    dump(ds, output)

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
def peptide_md_to_fasta(peptide_table, output_fasta):
    """
    convert peptide metadata to fasta format.

    For each peptide, we will add an entry to a fasta file with
    the unique peptide_id as the only entry into '>' header and
    oligonucleotide encoding on the line below.
    """
    convert_peptide_table_to_fasta(peptide_table, output_fasta)
