"""
@File: test_utils.py

@Author: Jared Galloway

Some unit tests for phippery.PhipData
"""

# built-in
import os
import sys

# dependency
import pytest
import numpy as np
import pandas as pd
import xarray as xr
import glob
from phippery.PhipData import load_counts
from phippery.PhipData import collect_sample_metadata
from phippery.PhipData import collect_peptide_metadata
from phippery.PhipData import collect_merge_prune_count_data


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


def test_load_counts(shared_datadir):

    pds = load_counts(
        counts_files=[
            (shared_datadir / f"test_files/counts/{tfile}")
            for tfile in os.listdir((shared_datadir / "test_files/counts"))
        ],
        peptide_metadata=(shared_datadir / "test_files/peptide_metadata.csv"),
        sample_metadata=(shared_datadir / "test_files/sample_metadata.csv"),
    )
    assert type(pds) == xr.Dataset
