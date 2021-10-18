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
from phippery.utils import id_coordinate_from_query
from phippery.utils import sample_id_coordinate_from_query
from phippery.utils import peptide_id_coordinate_from_query
from phippery.phipdata import convert_peptide_table_to_fasta
from phippery.phipdata import collect_merge_prune_count_data
from phippery.phipdata import add_stats
from phippery.phipdata import collect_sample_table
from phippery.phipdata import collect_peptide_table
from phippery.phipdata import dataset_from_csv
from phippery.phipdata import stitch_dataset
from phippery.phipdata import load
from phippery.phipdata import dump
from phippery.phipdata import dataset_to_wide_csv
# TODO move tidy to phipdata, yea?
from phippery.tidy import tidy_ds
from phippery.string import string_ds

from phippery.normalize import cpm


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
        click.echo(e)

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
        click.echo(e)

    # call backend for this one
    info = string_feature(
        ds,
        feature,
        dim=dimension
    )

    # write it
    click.echo(info)
    

@cli.command(name="query-expression")
@argument(
    'filename', 
    type=click.Path(exists=True)
)
@argument(
    'expression', 
    type=str
)
@click.option(
        '-d', '--dimension',
        type=click.Choice(['sample', 'peptide'], 
            case_sensitive=False),
        default='sample'
)
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
    )
def query(filename, expression, dimension, output):
    """
    Perform a single pandas-style query expression on dataset samples

    This command takes a single string query statement,
    aplied it to the sample table in the dataset,
    and returns the dataset with the slice applied.
    This mean that all enrichment layers get sliced.
    If no output (-o) is provided, by default this command
    will overwrite the provided dataset file. 

    Some rst:

    .. note:: for more information on pandas query style strings,
       please see the 
       `Pandas documentation 
       <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_    
       additionally, I've found `this blog
       <https://queirozf.com/entries/pandas-query-examples-sql-like-syntax-queries-in-dataframes>`_ 
       very helpful for performing queris on a dataframe. 
    """

    try:
        ds = load(filename)
    except Exception as e:
        click.echo(e)

    if output == None:
        if click.confirm(f'Without providing output path, you overwrite {filename}. Do you want to continue?'):
            output = filename
        else:
            click.echo('Abort')    
            return 
    
    if dimension == "sample":
        q = sample_id_coordinate_from_query(
            ds,
            query_list = [expression]
        )
    else:
        q = peptide_id_coordinate_from_query(
            ds,
            query_list = [expression]
        )

    dump(ds.loc[{f"{dimension}_id":q}], output)



@cli.command(name="query-table")
@argument('filename', type=click.Path(exists=True))
@argument('expression-table', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
    )
def query_table(filename, expression_table, output):
    """
    Perform dataset index a csv giving a set of query expressions

    This command takes a dsvpoviding a set
    of queries for both samples or peptide and
    applies each to the respective annotation table in the dataset,
    and returns the dataset with the slice applied.
    This mean that all enrichment layers get sliced.
    If no output (-o) is provided, by default this command
    will overwrite the provided dataset file. 

    An example csv might look like the following

    
    >dimension, expression
    >sample, "Cohort == 2.0"
    >sample, "technical_replicate_id > 500"

    Some rst:

    .. note:: for more information on pandas query style strings,
       please see the 
       `Pandas documentation 
       <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_
       additionally, I've found `this blog
       <https://queirozf.com/entries/pandas-query-examples-sql-like-syntax-queries-in-dataframes>`_ 
       very helpful for performing queris on a dataframe. 
    """

    try:
        ds = load(filename)
    except Exception as e:
        click.echo(e)
        return

    if output == None:
        if click.confirm(f'Without providing output path, you overwrite {filename}. Do you want to continue?'):
            output = filename
        else:
            click.echo('Abort')    
            return 

    query_df = pd.read_csv(expression_table, header=0)
    if len(query_df) == 0:
        click.echo("Error: empty expression table.")
        return

    click.echo(query_df.columns)
    for rc in ["dimension", "expression"]:
        if rc not in query_df.columns:
            click.echo(f"Error: csv must have {rc} column.")
            return

    sid, pid = id_coordinate_from_query(ds, query_df)
    if len(sid) == 0:
        click.echo("Error: query resulted in zero valid samples")
    if len(pid) == 0:
        click.echo("Error: query resulted in zero valid peptides")

    dump(ds.loc[{"sample_id":sid, "peptide_id":pid}], output)


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
    Export the given dataset to a tall style dataframe.
    """

    try:
        ds = load(filename)
    except Exception as e:
        click.echo(e)
        return

    tidy_ds(ds).to_csv(output)


# To wide
@cli.command(name="to-wide-csv")
@argument('filename', type=click.Path(exists=True))
@option(
    '-o','--output-prefix', 
    type=click.Path(exists=False), 
    default=None,
    required=True
    )
def to_wide_csv(filename, output_prefix):
    """
    Export the given dataset to wide style dataframes.
    """
    
    try:
        ds = load(filename)
    except Exception as e:
        click.echo(e)
        return
    
    dataset_to_wide_csv(ds, output_prefix)




# Fit Gam-Pois
@cli.command(name="fit-predict-neg-binom")
@argument('fitting-dataset', type=click.Path(exists=True))
@argument('predict-dataset', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
)
def fit_predict_neg_binom(fitting_dataset, predict_dataset, output):
    """
    STUB - coming soon

    Fit and predict mlxp values using the negative binomial model
    """

    try:
        ds = load(filename)
    except Exception as e:
        click.echo(e)
    pass

# Calculate fold enrichment
@cli.command(name="fold-enrichment")
@argument('dataset', type=click.Path(exists=True))
@argument('control-dataset', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
)
def fold_enrichment(fitting_dataset, conrol_dataset, output):
    """
    STUD - coming soon
    """


# Calculate fold enrichment
@cli.command(name="cpm")
@argument('dataset', type=click.Path(exists=True))
@option(
    '-o','--output', 
    type=click.Path(exists=False), 
    default=None,
    required=False
)
def cpm(fitting_dataset, conrol_dataset, output):
    """
    STUB - coming soon
    """
    
## Calculate fold enrichment
#@cli.command(name="diff")
#@argument('dataset', type=click.Path(exists=True))
#@option(
#    '-o','--output', 
#    type=click.Path(exists=False), 
#    default=None,
#    required=False
#)
#def cpm(fitting_dataset, conrol_dataset, output):
#    """
#    STUB - coming soon
#    """
#    
## Calculate fold enrichment
#@cli.command(name="cpm")
#@argument('dataset', type=click.Path(exists=True))
#@option(
#    '-o','--output', 
#    type=click.Path(exists=False), 
#    default=None,
#    required=False
#)
#def cpm(fitting_dataset, conrol_dataset, output):
#    """
#    STUB - coming soon
#    """
#    

# TODO Enrichment
# TODO Diff Sel 

