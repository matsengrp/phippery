"""
@File: test_phip_data.py

@Author: Jared Galloway

A set of unit tests for the gathering and exporting of phip data.
"""

# built-in dependency
import os
import sys
import pickle
import glob

# dependency
import pytest
import numpy as np
import pandas as pd
import xarray as xr

from sim_test_generator import SimulationTest
from sim_test_generator import iter_sim_tests
from sim_test_generator import generate_sim_ds


# local dependency
from phippery.phipdata import counts_metadata_to_dataset
from phippery.phipdata import collect_sample_metadata
from phippery.phipdata import collect_peptide_metadata
from phippery.phipdata import dataset_to_csv
from phippery.phipdata import csv_to_dataset
from phippery.phipdata import df_to_dataset
from phippery.phipdata import load
from phippery.phipdata import dump
from phippery.phipdata import add_stats


def test_add_stats(shared_datadir, tmp_path):
    """

    """
    ds = generate_sim_ds()
    d = tmp_path / "sub"
    d.mkdir()
    files = []
    for sid in ds.sample_id.values:
        fp = open(f"{d}/{sid}.txt", "w")
        files.append(f"{d}/{sid}.txt")
        for stat in ["stat_a", "stat_b", "stat_c"]:
            fp.write(f"{stat}\t{np.random.randint(10)}\n")
        fp.close()
    ds = add_stats(ds, files)


# TODO
def test_load(shared_datadir):
    """
    simple wrapper for loading xarray datasets from pickle binary
    """

    pass


# TODO
def test_dump(shared_datadir):
    """
    simple wrapper for dump'ing xarray datasets to pickle binary
    """

    pass


def test_sims_generator(shared_datadir):

    for sim_test in iter_sim_tests(shared_datadir):
        assert type(sim_test) == SimulationTest
        assert type(sim_test.counts) == list


def test_counts_metadata_to_dataset(shared_datadir):

    for sim_test in iter_sim_tests(shared_datadir):
        pds = counts_metadata_to_dataset(
            sim_test.counts, sim_test.pep_meta, sim_test.sam_meta
        )
        assert type(pds) == xr.Dataset
        assert np.all(pds.counts.values == sim_test.solution)
        # assert pds.attrs["sample_coord_dim"] == "sample_id"
        # assert pds.attrs["peptide_coord_dim"] == "peptide_id"


def test_read_write_csv(shared_datadir, tmp_path):

    for sim_test in iter_sim_tests(shared_datadir):
        d = tmp_path / "sub"
        d.mkdir()
        # p = d / "test_"
        dataset_to_csv(sim_test.pds, f"{d}/test")
        assert os.path.exists(f"{d}/test_counts.csv")
        assert os.path.exists(f"{d}/test_sample_table.csv")
        assert os.path.exists(d / "test_peptide_table.csv")

        ds_from_csv = csv_to_dataset(
            f"{d}/test_counts.csv",
            f"{d}/test_peptide_table.csv",
            f"{d}/test_sample_table.csv",
        )

        assert type(ds_from_csv) == xr.Dataset
        assert sim_test.pds.equals(ds_from_csv)


def test_df_to_dataset(shared_datadir):

    for sim_test in iter_sim_tests(shared_datadir):
        ds = df_to_dataset(
            counts_df=sim_test.pds.counts.to_pandas(),
            peptide_table_df=sim_test.pds.peptide_table.to_pandas(),
            sample_table_df=sim_test.pds.sample_table.to_pandas(),
        )
        assert sim_test.pds.equals(ds)
