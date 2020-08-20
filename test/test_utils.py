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
from phippery.utils import get_all_sample_metadata_factors
from phippery.utils import get_all_peptide_metadata_factors


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


def get_test_pds(shared_datadir):
    """
    a helper function to collect counts into an xarray dataset for testing.
    """

    pds = counts_metadata_to_dataset(
        counts_files=[
            (shared_datadir / f"test_files/counts/{tfile}")
            for tfile in os.listdir((shared_datadir / "test_files/counts"))
        ],
        peptide_metadata=(shared_datadir / "test_files/peptide_metadata.csv"),
        sample_metadata=(shared_datadir / "test_files/sample_metadata.csv"),
    )

    return pds


def test_get_all_metadata_factors(shared_datadir):
    """
    Test the get_all_x_factors utilities.
    """

    pds = get_test_pds(shared_datadir)
    assert len(get_all_sample_metadata_factors(pds, "reference")) == 1
    assert len(get_all_peptide_metadata_factors(pds, "Oligo")) == 10
