"""
@File: phippery_cli.py

@Author: Jared Galloway

Command line interface (CLI) for phippery.
"""

# built-in
import gzip
import glob

# dependencies
import pandas as pd
import xarray as xr
from click import Path, group, option, argument
import click

# local
from phippery import utils
from phippery.string import string_ds


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
    sample_table,
    peptide_table,
    counts_matrix,
    output,
):
    """
    Load and dump xarray dataset given a set of wide csv's

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

    ds = utils.dataset_from_csv(counts_matrix, peptide_table, sample_table)
    utils.dump(ds, output)


@cli.command(name="about")
@argument("filename", type=click.Path(exists=True))
@option("-v", "--verbose", count=True)
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

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)

    # call backend for this one
    info = string_ds(ds, verbose)

    # write it
    click.echo(info)


from phippery.string import string_feature


@cli.command(name="about-feature")
@argument("feature", type=str)
@argument("filename", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dimension",
    type=click.Choice(["sample", "peptide"], case_sensitive=False),
    default="sample",
    help="The dimension we expect to find this feature"
)
@click.option(
    "--distribution/--counts",
    default=True,
    help="Force a specific output of either value counts or distribution for quantitative features"
)
# def string_feature(ds, feature: str, verbosity = 0, dim="sample"):
def about_feature(filename, dimension, feature, distribution):
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
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)

    # call backend for this one
    info = string_feature(ds, feature, dim=dimension, numeric_dis=distribution)

    # write it
    click.echo(info)


@cli.command(name="split-groups")
@click.option(
    "-d",
    "--dimension",
    type=click.Choice(["sample", "peptide"], case_sensitive=False),
    default="sample",
)
@option("split-features", type=str)
@argument("filename", type=click.Path(exists=True))
def query_expression(filename, expression, dimension, output):
    """
    Perform a single pandas-style query expression on dataset samples

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
       additionally, I've found `this blog
       <https://queirozf.com/entries/pandas-query-examples-sql-like-syntax-queries-in-dataframes>`_
       very helpful for performing queries on a dataframe.
    """

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)

    if output == None:
        output = "sliced_dataset.phip"

    q = utils.id_query(ds, expression, dimension)

    utils.dump(ds.loc[{f"{dimension}_id": q}], output)


@cli.command(name="merge")
@option("-o", "--output", type=click.Path(exists=False), default=None, required=False)
@argument(
    "datasets",
)
def merge(output, datasets):
    """ """

    try:
        dss = [utils.load(f) for f in glob.glob(datasets)]
    except Exception as e:
        click.echo(e)

    if output == None:
        output = "merged_dataset.phip"

    merged = xr.merge(dss)
    utils.dump(merged, output)


##################
##################


@cli.command(name="query-expression")
@click.option(
    "-d",
    "--dimension",
    type=click.Choice(["sample", "peptide"], case_sensitive=False),
    default="sample",
)
@option("-o", "--output", type=click.Path(exists=False), default=None, required=False)
@argument("expression", type=str)
@argument("filename", type=click.Path(exists=True))
def query_expression(filename, expression, dimension, output):
    """
    Perform a single pandas-style query expression on dataset samples

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
       additionally, I've found `this blog
       <https://queirozf.com/entries/pandas-query-examples-sql-like-syntax-queries-in-dataframes>`_
       very helpful for performing queries on a dataframe.
    """

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)

    if output == None:
        output = "sliced_dataset.phip"

    q = utils.id_query(ds, expression, dimension)

    utils.dump(ds.loc[{f"{dimension}_id": q}], output)


@cli.command(name="query-table")
@argument("filename", type=click.Path(exists=True))
@argument("expression-table", type=click.Path(exists=True))
@option("-o", "--output", type=click.Path(exists=False), default=None, required=False)
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

    \f

    An example table (csv) might look like the following

    .. list-table:: Example
        :widths: 25 25
        :header-rows: 1

        * - dimension 
          - expression
        * - sample
          - "Cohort == 2.0"
        * - sample
          - "technical_replicate_id > 500"

    .. note:: for more information on pandas query style strings,
       please see the
       `Pandas documentation
       <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_
       additionally, I've found `this blog
       <https://queirozf.com/entries/pandas-query-examples-sql-like-syntax-queries-in-dataframes>`_
       very helpful for performing queries on a dataframe.
    """

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)
        return

    if output == None:
        if click.confirm(
            f"Without providing output path, you overwrite {filename}. Do you want to continue?"
        ):
            output = filename
        else:
            click.echo("Abort")
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

    utils.dump(ds.loc[{"sample_id": sid, "peptide_id": pid}], output)


# To tall
@cli.command(name="to-tall-csv")
@argument("filename", type=click.Path(exists=True))
@option("-o", "--output", type=click.Path(exists=False), default=None, required=True)
def to_tall_csv(filename: str, output: str):
    """
    Export the given dataset to a tall style dataframe.
    """

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)
        return

    # Open a connection to the output file
    if output.endswith(".gz"):
        handle = gzip.open(output, 'wt')
    else:
        handle = open(output, 'w')

    # Generate tall tables for each of the samples in turn
    # Each of the tables for each sample will have the same column order
    for i, sample_df in enumerate(utils.yield_tall(ds)):

        # If it's the first one, include the header
        handle.write(
            sample_df.to_csv(
                header=i == 0,
                index=False,
                na_rep="NA"
            )
        )
    handle.close()


# To wide
@cli.command(name="to-wide-csv")
@argument("filename", type=click.Path(exists=True))
@option(
    "-o", "--output-prefix", type=click.Path(exists=False), default=None, required=True
)
def to_wide_csv(filename, output_prefix):
    """
    Export the given dataset to wide style dataframes.
    """

    try:
        ds = utils.load(filename)
    except Exception as e:
        click.echo(e)
        return

    utils.to_wide_csv(ds, output_prefix)
