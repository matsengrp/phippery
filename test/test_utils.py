"""
@File: test_utils.py

@Author: Jared Galloway

Some unit tests for phippery.utils
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
import phippery.utils as utils


def test_files(shared_datadir):
    """
    simple test to make sure we can find our test files,

    (shared_datadir/"test_files"/"counts") - is a directory,
    with  7 files:

    1.1.tsv  1.2.tsv  2.1.tsv  2.2.tsv  35.1.tsv  35.2.tsv  37.1.tsv

    each a name for a file containing example peptide alignment counts for
    a sample in a phip-seq experiement. The first integer of the filename
    is the sample number and the second number is th technical replicate.

    Each file contains two columns seperated by tabs. first column
    is peptide id's and the second is the alignment counts.
    """

    # assert os.path.exists(self.counts)
    # assert os.path.exists(self.sample_metadata)
    # assert os.path.exists(self.peptide_metadata)

    assert os.path.exists((shared_datadir / "test_files/counts"))
    assert os.path.exists((shared_datadir / "test_files" / "sample_metadata.tsv"))
    assert os.path.exists((shared_datadir / "test_files" / "peptide_metadata.tsv"))

    return None


def test_collect_sample_metadata_types(shared_datadir):
    """test that we're getting the right types"""

    sample_md = utils.collect_sample_metadata(
        (shared_datadir / "test_files/sample_metadata.tsv")
    )
    assert type(sample_md) == xr.DataArray


def test_collect_peptide_metadata_types(shared_datadir):
    """test that we're getting the right types"""

    peptide_md = utils.collect_peptide_metadata(
        (shared_datadir / "test_files/peptide_metadata.tsv")
    )
    assert type(peptide_md) == xr.DataArray


# def test_shared_dict(shared_datadir):
#    contents = (shared_datadir / "test_files/meh.txt").read_text()
#    assert contents == 'ey\n'


def test_collect_merge_prune_count_data(shared_datadir):
    """test that we're getting the right types"""

    # print((shared_datadir / "test_files/counts"))
    # print(glob.glob("test/data/test_files/counts/*.tsv"))
    counts = utils.collect_merge_prune_count_data(
        (shared_datadir / "test_files/counts")
    )

    assert type(counts) == xr.DataArray
