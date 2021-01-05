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

# functions to test
from phippery.utils import trim_index
from phippery.utils import convert_peptide_metadata_to_fasta
from phippery.utils import get_all_sample_metadata_factors
from phippery.utils import get_all_peptide_metadata_factors
from phippery.utils import iter_sample_groups
from phippery.utils import iter_peptide_groups
from phippery.utils import id_coordinate_subset


def test_trim_index():
    """
    Test trimming function
    """

    assert trim_index("aaaTTTaaa") == "TTT"
    assert trim_index("") == ""
    assert trim_index("TTT") == "TTT"
    assert trim_index("aaaTTT") == "TTT"
    assert trim_index("TTTaaa") == "TTT"


# TODO simulate a dataset directly to test these, it'll be a lot faster if we keep things
# in memory rather than reading in the simulation tests each time.
def test_convert_peptide_metadata_to_fasta(shared_datadir, tmp_path):
    """
    Test convert to fasta from peptide metadata.

    Given a unique identifier and a set of Oligos, we expect the fasta file
    to have a header ">" and Oligo sequence for each.
    """

    for sim_test in iter_sim_tests(shared_datadir):

        d = tmp_path / "conv"
        d.mkdir()
        p = d / "test_fasta"
        convert_peptide_metadata_to_fasta(sim_test.pep_meta, p)
        assert os.path.exists(p)
        pep_meta = sim_test.pds.peptide_table
        for i, line in enumerate(open(p, "r")):
            line = line.strip()
            if i % 2 == 0:
                assert line.startswith(">")
                assert int(line[1:]) == pep_meta.loc[i // 2, :].peptide_id
            else:
                pep_meta_oligo = trim_index(str(pep_meta.loc[i // 2, "Oligo"].values))
                assert line == pep_meta_oligo


def test_get_all_sample_metadata_factors(shared_datadir):
    """
    test to get all factors from sample table feature
    """

    for sim_test in iter_sim_tests(shared_datadir):

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

                # assert the xarray dataset subset has the same features and strucure
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
                # assert we hit all groups
                assert group in peptide_table.loc[:, feature].values

                # assert the xarray dataset subset has the same features and strucure
                assert np.all(
                    group_ds.data_vars.keys() == sim_test.pds.data_vars.keys()
                )
                assert np.all(group_ds.attrs.keys() == sim_test.pds.attrs.keys())
                assert np.all(group_ds.coords.keys() == sim_test.pds.coords.keys())


def test_id_coordinate_subset(shared_datadir):
    """
    test getting sample coordinates from subset function
    by comparing the indexing to pandas series indexing
    """

    for sim_test in iter_sim_tests(shared_datadir):
        for feat in sim_test.pds.sample_metadata.values:
            possible_fact = get_all_sample_metadata_factors(sim_test.pds, feat)
            for fac in possible_fact:

                sample_series = sim_test.pds.sample_table.to_pandas().loc[:, feat]

                # testing is_equal_to
                sam_subset = id_coordinate_subset(
                    sim_test.pds, where=feat, is_equal_to=fac
                )

                sam_subset_series = list(sample_series[sample_series == fac].index)
                assert np.all(sam_subset == sam_subset_series)

                # testing not equal to
                sam_subset = id_coordinate_subset(
                    sim_test.pds, where=feat, is_not_equal_to=fac
                )
                sam_subset_series = list(sample_series[sample_series != fac].index)
                assert np.all(sam_subset == sam_subset_series)

                # testing is_in
                sam_subset = id_coordinate_subset(
                    sim_test.pds, where=feat, is_in=possible_fact
                )
                assert len(sam_subset) == len(sim_test.pds.sample_id.values)


def test_id_coordinate_subset_table_error(shared_datadir):
    for sim_test in iter_sim_tests(shared_datadir):
        with pytest.raises(ValueError):
            id_coordinate_subset(sim_test.pds, where="foo", table="boo")


def test_id_coordinate_subset_table_error_table(shared_datadir):
    for sim_test in iter_sim_tests(shared_datadir):
        with pytest.raises(ValueError):
            id_coordinate_subset(sim_test.pds, where="foo")


def test_id_coordinate_subset_table_error_where(shared_datadir):
    for sim_test in iter_sim_tests(shared_datadir):
        with pytest.raises(ValueError):
            id_coordinate_subset(
                sim_test.pds, where=sim_test.pds.sample_metadata.values[0]
            )
