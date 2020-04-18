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
import xarray as xr
import phippery.utils as utils


class TestUtils:
    """ Object to test phippery.utils """

    # TODO
    # can shared_datadir be accessed without passing
    # it as a parameter?
    # self.counts = (shared_datadir/"test_files"/"counts")
    # self.sample_metadata = (shared_datadir/"test_files"/"sample_metadata.tsv")
    # self.peptide_metadata = (shared_datadir/"test_files"/"peptide_metadata.tsv")

    def test_files(self, shared_datadir):
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

        assert os.path.exists((shared_datadir / "test_files" / "counts"))
        assert os.path.exists((shared_datadir / "test_files" / "sample_metadata.tsv"))
        assert os.path.exists((shared_datadir / "test_files" / "peptide_metadata.tsv"))

        return None

    def test_collect_sample_metadata_types(self, shared_datadir):
        """test that we're getting the right types"""

        # sample_md = utils.collect_sample_metadata(
        #    f"{shared_datadir}/test_files/sample_metadata.tsv"
        # )

        # assert type(sample_md) == xr.DataArray
