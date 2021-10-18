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
from click import Choice, Path, command, group, option, argument
import click

# local
from phippery.phipdata import convert_peptide_table_to_fasta
from phippery.phipdata import collect_merge_prune_count_data
from phippery.phipdata import add_stats
from phippery.phipdata import collect_sample_table
from phippery.phipdata import collect_peptide_table
from phippery.phipdata import dataset_from_csv
from phippery.phipdata import stitch_dataset
from phippery.phipdata import load
from phippery.phipdata import dump
from phippery.string import string_ds


# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def phipflowcli():
    """
    A set of tools for aiding in the pre-processing
    alignment pipeline.

    These are NOT recommended to be used in downstream analysis
    without extreme caution.
    """
    pass

@phipflowcli.command(name="load-from-counts-tsv")
@option(
    "-s",
    "--sample_table",
    required=True,
    type=Path(exists=True),
    help="Path to sample table csv.",
)
@option(
    "-p",
    "--peptide_table",
    required=True,
    type=Path(exists=True),
    help="Path to peptide table csv.",
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
    required=False,
    help="File pattern (glob) for all counts files to be merged",
)
@option(
    "-o",
    "--output",
    required=True,
    help="Path where the phip dataset will be dump'd to netCDF",
)
def load_from_counts_tsv(
    sample_table, 
    peptide_table, 
    counts_file_pattern, 
    stats_file_pattern, 
    output
):
    """
    (Primarily) PhIP-Flow pipeline helper function

    Collect sample and peptide metadata tables along with a
    two-column tsv file for each sample, and produce a properly formatted xarray dataset.
    """

    counts = [f for f in glob.glob(counts_file_pattern)]
    stats = [f for f in glob.glob(stats_file_pattern)]

    merged_counts = collect_merge_prune_count_data(counts)
    peptide_table = collect_peptide_table(peptide_table)
    sample_table = collect_sample_table(sample_table)

    ds = stitch_dataset(
        counts=merged_counts, peptide_table=peptide_table, sample_table=sample_table,
    )

    # TODO Finish add stats check 
    # ds = add_stats(ds, stats)
    dump(ds, output)


@phipflowcli.command(name="peptide-md-to-fasta")
@option(
    "-d",
    "--peptide-table",
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
    (Primarily) PhIP-Flow helper function

    convert peptide metadata to fasta format.
    For each peptide, we will add an entry to a fasta file with
    the unique peptide_id as the only entry into '>' header and
    oligonucleotide encoding on the line below.
    """
    convert_peptide_table_to_fasta(peptide_table, output_fasta)



# entry point
@group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Welcome to the phippery CLI! 

    Here we present a few commands that allow users to
    slice, transform, normalize, fit models to, and more given 
    a binary pickle dump'd xarray, usually as a result of running the
    PhIP-Flow pipeline. 

    For more information and example workflows please refer to 
    the full documentation 
    at https://matsengrp.github.io/phippery/
    """
    pass


# TODO the sample and peptide tables aren't required.
@cli.command(name="load-from-csv")
@option(
    "-s",
    "--sample_table",
    required=True,
    type=Path(exists=True),
    help="Path to sample table csv.",
)
@option(
    "-p",
    "--peptide_table",
    required=True,
    type=Path(exists=True),
    help="Path to peptide table csv.",
)
@option(
    "-c",
    "--counts_matrix",
    required=True,
    type=Path(exists=True),
    help="Path to counts matrix csv.",
)
@option(
    "-o",
    "--output",
    required=True,
    help="Path where the phip dataset will be dump'd to netCDF",
)
def load_from_csv(
    sample_table, peptide_table, counts_matrix, output,
):
    """
    Load and create pickle dump'd xarray dataset from a set of wide csv's 
    assuming they have the correct shapes and format. 
   
    Using this command usually means you have either:
    
    1. Decided to store the output of your analysis in the form
       of wide csv's instead of a pickle dump'd binary for 
       longer-term storage.

    2. Created your own enrichment matrix
       without the help of the phip-flow alignment pipeline.
    \f

    .. note::
      In the case of #2, please note that your matrix data
      must be numeric and have shape (len(peptide_table), len(sample_table)).
      Finally, you must include
      :ref:`pipeline outputs <sec_pipeline_outputs>`

    .. note::
      Currently only accepting a single enrichment matrix.
    """

    ds = dataset_from_csv(
        counts_matrix,
        peptide_table,
        sample_table
    )
    dump(ds, output)


@cli.command(name="about")
@argument('filename', type=click.Path(exists=True))
@option('-v', '--verbose', count=True)
def about(filename, verbose):
    """
    Summarize the data in a given dataset

    If no verbosity flag is set, this will print the 
    basic information about number of
    samples, peptides, and enrichment layers
    in a given dataset. With a verbosity of one (-v) you will
    get also information about annotations and available datatypes.
    If verbosity flag is set to two (-vv) - Print 
    detailed information about all data tables including annotation
    feature types and distribution of enrichment values for each layer.
    A verbosity of three will basically loop through all features
    and give you you the detailed description of each.
    """

    # TODO Here I should be calling some sort of validation
    # probably should have these checks in load()
    try:
        ds = load(filename)
    except Exception as e:
        st.write(e)

    # call backend for this one
    info = string_ds(ds, verbose)
    
    # write it
    click.echo(info)


from phippery.string import string_feature
# TODO verbosity
@cli.command(name="about-feature")
@argument('filename', type=click.Path(exists=True))
@click.option(
        '-d', '--dimension',
        type=click.Choice(['sample', 'peptide'], 
            case_sensitive=False),
        default='sample'
)
@argument('feature', type=str)
#def string_feature(ds, feature: str, verbosity = 0, dim="sample"):
def about_feature(filename, dimension, feature):
    """
    Summarize details about a specific sample or peptide annotation feature.

    The function will tell you information about a specific feature
    in you sample annnotation table, depending on it's inferred datatype.
    For numeric feature types the command will get information about quantiles,
    for categorical or boolean feature types, the function will give
    individual factor-level counts.

    Both datatype categories will print a few examples of valid dataset
    queries using the feature in question
    """

    try:
        ds = load(filename)
    except Exception as e:
        st.write(e)

    # call backend for this one
    info = string_feature(
        ds,
        feature,
        dim=dimension
    )

    # write it
    click.echo(info)
    


#@cli.command(name="about-pepide-feature")
#@argument('filename', type=click.Path(exists=True))
#@argument('feature', type=str)
#def about_peptide_feature(filename, feature):
#    """
#    STUB - Summarize details about a specific peptide annotation feature.
#
#    The function will infer telll you information about a specific feature
#    in you peptide annnotation table, depending on it's inferred datatype.
#    For numeric feature types the command will get information about quantiles,
#    for categorical or boolean feature types, the function will give
#    individual factor-level counts.
#
#    Both datatype categories will print a few examples of valid dataset
#    queries using the feature in question
#    """
#    try:
#        ds = load(filename)
#    except Exception as e:
#        st.write(e)
#    # call backend for this one
#    #info = about_phipdata_peptide_feature(ds, verbosity)
#
#    # write it
#    #st.write(info)
#    pass


@cli.command(name="query-samples")
@argument(
    'filename', 
    type=click.Path(exists=True)
)
@argument(
    'expression', 
    type=str
)
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
    )
def query_samples(filename, expression):
    """
    STUB - Perform a single pandas-style query expression on dataset samples

    This command takes a single string query statement,
    aplied it to the sample table in the dataset,
    and returns the dataset with the slice applied.
    This mean that all enrichment layers get sliced.
    If no output (-o) is provided, by default this command
    will overwrite the provided dataset file. 
    \f

    .. note:: for more information on pandas query style strings,
       please see the 
       `Pandas documentation 
       <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_

    
    """
    pass


@cli.command(name="query-table")
@argument('filename', type=click.Path(exists=True))
@argument('expression-table', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
    )
def query_table(filename, expression_table):
    """
    STUB - Perform a single pandas-style query expression on dataset table

    This command takes a single string query statement,
    aplied it to the sample table in the dataset,
    and returns the dataset with the slice applied.
    This mean that all enrichment layers get sliced.
    If no output (-o) is provided, by default this command
    will overwrite the provided dataset file. 
    \f

    .. note:: for more information on pandas query style strings,
       please see the 
       `Pandas documentation 
       <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_

    
    """
    pass

# To tall
@cli.command(name="to-tall-csv")
@argument('filename', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=True
    )
def to_tall_csv(filename, output):
    """
    STUB - Export the given dataset to a Tall style dataframe.
    """
    pass


# To wide
@cli.command(name="to-wide-csv")
@argument('filename', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=True
    )
def to_wide_csv(filename, output):
    """
    STUB - Export the given dataset to a Tall style dataframe.
    """
    pass


# Fit Gam-Pois
@cli.command(name="fit-neg-binom")
@argument('fitting-dataset', type=click.Path(exists=True))
@argument('predict-dataset', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=True
    )
def fit_neg_binom(fitting_dataset, predict_dataset, output):
    """
    STUB - Fit and predict mlxp values using the negative binomial model
    """
    pass


# Enrichment



