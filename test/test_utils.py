"""
@File: test_utils.py

@Author: Jared Galloway

Some unit tests for phippery.utils
"""

# dependency
import pytest
import phippery.utils


def load_test_files():
    """
    This function simply gives the paths
    to the following three test files.

    counts: This is a directory with some
        tab seperated value alignment counts
        for demultiplexed samples:
        1.1, 1.2, 2.1, 2.2, 35.1, 35.2, 37.1
        each with a count for peptides:
        1,2,3,4,5.

    sample_metadata:
    """
    test_files = "../data/test_files/"
    return (
        f"{test_files}counts/",
        f"{test_files}sample_metadata",
        f"{test_files}peptide_metadata",
    )


def test_collect_merge_prune_data():
    """
    simple test to get things up and running
    """
    # counts, s_meta, p_meta =
    assert True
