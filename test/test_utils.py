"""
@File: test_phip_data.py

@Author: Jared Galloway

A set of unit tests for the gathering and exporting of phip data.
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

from sim_test_generator import SimulationTest
from sim_test_generator import iter_sim_tests
from sim_test_generator import make_hardcoded_ds

# functions to test
from phippery.utils import get_all_sample_metadata_factors
from phippery.utils import get_all_peptide_metadata_factors
from phippery.utils import iter_sample_groups
from phippery.utils import iter_peptide_groups


def test_get_all_sample_metadata_factors(shared_datadir):
    """
    test to get all factors from sample table feature
    """

    for sim_test in iter_sim_tests(shared_datadir):
        print(sim_test)

        sample_table = sim_test.pds.sample_table.to_pandas().reset_index()
        for feature in sim_test.pds.sample_metadata.values:
            factors = get_all_sample_metadata_factors(sim_test.pds, feature)
            assert np.all([f in sample_table.loc[:, feature]] for f in factors)


def test_get_all_peptide_metadata_factors(shared_datadir):
    """
    test to get all factors from peptide table feature
    """

    for sim_test in iter_sim_tests(shared_datadir):

        peptide_table = sim_test.pds.peptide_table.to_pandas().reset_index()
        for feature in sim_test.pds.peptide_metadata.values:
            factors = get_all_peptide_metadata_factors(sim_test.pds, feature)
            assert np.all([f in peptide_table.loc[:, feature]] for f in factors)


def test_iter_sample_groups(shared_datadir):
    """
    test iteration of sample groups
    """

    for sim_test in iter_sim_tests(shared_datadir):

        sample_table = sim_test.pds.sample_table.to_pandas().reset_index()
        for feature in sim_test.pds.sample_metadata.values:
            for group, group_ds in iter_sample_groups(sim_test.pds, feature):
                # assert we hit all groups
                assert group in sample_table.loc[:, feature].values
                assert np.all(
                    group_ds.data_vars.keys() == sim_test.pds.data_vars.keys()
                )
                assert np.all(group_ds.attrs.keys() == sim_test.pds.attrs.keys())
                assert np.all(group_ds.coords.keys() == sim_test.pds.coords.keys())


def test_iter_peptide_groups(shared_datadir):
    """
    test iteration of peptide groups
    """

    for sim_test in iter_sim_tests(shared_datadir):

        peptide_table = sim_test.pds.peptide_table.to_pandas().reset_index()
        for feature in sim_test.pds.peptide_metadata.values:
            for group, group_ds in iter_peptide_groups(sim_test.pds, feature):
                assert group in peptide_table.loc[:, feature].values
                assert np.all(
                    group_ds.data_vars.keys() == sim_test.pds.data_vars.keys()
                )
                assert np.all(group_ds.attrs.keys() == sim_test.pds.attrs.keys())
                assert np.all(group_ds.coords.keys() == sim_test.pds.coords.keys())
