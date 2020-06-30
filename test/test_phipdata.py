"""
@File: test_utils.py

@Author: Jared Galloway

Some unit tests for phippery.phipdata
"""

# built-in dependency
import os
import sys

# dependency
import pytest
import numpy as np
import pandas as pd
import xarray as xr
import glob

# local dependency
from phippery.phipdata import counts_metadata_to_dataset
from phippery.phipdata import collect_sample_metadata
from phippery.phipdata import collect_peptide_metadata
from phippery.phipdata import collect_merge_prune_count_data
from phippery.phipdata import dataset_to_csv
from phippery.phipdata import csv_to_dataset


def test_files(shared_datadir):
    """
    simple test to make sure we can find our test files,

    each a name for a file containing example peptide alignment counts for
    a sample in a phip-seq experiement. The first integer of the filename
    is the sample number and the second number is th technical replicate.

    Each file contains two columns seperated by tabs. first column
    is peptide id's and the second is the alignment counts.
    """

    assert os.path.exists((shared_datadir / "test_files/counts"))
    assert os.path.exists((shared_datadir / "test_files" / "sample_metadata.csv"))
    assert os.path.exists((shared_datadir / "test_files" / "peptide_metadata.csv"))


def test_collect_sample_metadata_types(shared_datadir):
    """test that we're getting the right types"""

    sample_md = collect_sample_metadata(
        (shared_datadir / "test_files/sample_metadata.csv")
    )
    assert type(sample_md) == pd.DataFrame


def test_collect_peptide_metadata_types(shared_datadir):
    """test that we're getting the right types"""

    peptide_md = collect_peptide_metadata(
        (shared_datadir / "test_files/peptide_metadata.csv")
    )
    assert type(peptide_md) == pd.DataFrame


def test_collect_merge_prune_count_data(shared_datadir):
    """test that we're getting the right types"""

    counts = collect_merge_prune_count_data(
        [
            (shared_datadir / f"test_files/counts/{tfile}")
            for tfile in os.listdir((shared_datadir / "test_files/counts"))
        ]
    )
    # first is counts df
    assert type(counts[0]) == pd.DataFrame
    # second is replicate info
    assert type(counts[1]) == pd.DataFrame


def test_counts_metadata_to_dataset(shared_datadir):

    pds = counts_metadata_to_dataset(
        counts_files=[
            (shared_datadir / f"test_files/counts/{tfile}")
            for tfile in os.listdir((shared_datadir / "test_files/counts"))
        ],
        peptide_metadata=(shared_datadir / "test_files/peptide_metadata.csv"),
        sample_metadata=(shared_datadir / "test_files/sample_metadata.csv"),
    )
    assert type(pds) == xr.Dataset


def test_dump_load_csv(shared_datadir, tmp_path):

    pds = counts_metadata_to_dataset(
        counts_files=[
            (shared_datadir / f"test_files/counts/{tfile}")
            for tfile in os.listdir((shared_datadir / "test_files/counts"))
        ],
        peptide_metadata=(shared_datadir / "test_files/peptide_metadata.csv"),
        sample_metadata=(shared_datadir / "test_files/sample_metadata.csv"),
    )

    d = tmp_path / "csv"
    d.mkdir()
    p = d / "test"

    dataset_to_csv(pds, p)
    pds_from_csv = csv_to_dataset(
        counts=f"{p}_counts.csv",
        peptide_table=f"{p}_peptide_table.csv",
        sample_table=f"{p}_sample_table.csv",
    )

    assert pds.equals(pds_from_csv)
