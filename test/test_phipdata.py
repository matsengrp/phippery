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
from phippery.phipdata import convert_peptide_table_to_fasta
from phippery.phipdata import stitch_dataset
from phippery.phipdata import collect_sample_table
from phippery.phipdata import collect_peptide_table
from phippery.phipdata import dataset_to_csv
from phippery.phipdata import load
from phippery.phipdata import dump
from phippery.phipdata import add_stats
from phippery.phipdata import trim_index


# TODO simulate a dataset directly to test these, it'll be a lot faster if we keep things
# in memory rather than reading in the simulation tests each time.
def test_convert_peptide_table_to_fasta(shared_datadir, tmp_path):
    """
    Test convert to fasta from peptide metadata.

    Given a unique identifier and a set of Oligos, we expect the fasta file
    to have a header ">" and Oligo sequence for each.
    """

    for sim_test in iter_sim_tests(shared_datadir):

        d = tmp_path / "conv"
        d.mkdir()
        p = d / "test_fasta"
        convert_peptide_table_to_fasta(sim_test.pep_meta, p)
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


def test_trim_index():
    """
    Test trimming function
    """

    assert trim_index("aaaTTTaaa") == "TTT"
    assert trim_index("") == ""
    assert trim_index("TTT") == "TTT"
    assert trim_index("aaaTTT") == "TTT"
    assert trim_index("TTTaaa") == "TTT"


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


def test_stitch_dataset(shared_datadir):

    for sim_test in iter_sim_tests(shared_datadir):

        pds = stitch_dataset(
            sim_test.pds.counts.to_pandas(),
            sim_test.pds.peptide_table.to_pandas(),
            sim_test.pds.sample_table.to_pandas(),
        )
        assert type(pds) == xr.Dataset
        assert np.all(pds.counts.values == sim_test.solution)


def test_read_write_csv(shared_datadir, tmp_path):

    load = lambda p: pd.read_csv(p, index_col=0)  # noqa
    for sim_test in iter_sim_tests(shared_datadir):
        d = tmp_path / "sub"
        d.mkdir()
        dataset_to_csv(sim_test.pds, f"{d}/test")
        assert os.path.exists(f"{d}/test_counts.csv")
        assert os.path.exists(f"{d}/test_sample_table.csv")
        assert os.path.exists(d / "test_peptide_table.csv")

        # print(load(f"{d}/test_counts.csv"))
        # print(load(f"{d}/test_peptide_table.csv"))
        # print(load(f"{d}/test_sample_table.csv"))

        counts = load(f"{d}/test_counts.csv")
        counts.index = counts.index.astype(int)
        counts.columns = counts.columns.astype(int)

        ds_from_csv = stitch_dataset(
            counts,
            collect_peptide_table(f"{d}/test_peptide_table.csv"),
            collect_sample_table(f"{d}/test_sample_table.csv"),
        )

        assert type(ds_from_csv) == xr.Dataset
        assert sim_test.pds.equals(ds_from_csv)


def test_df_to_dataset(shared_datadir):

    for sim_test in iter_sim_tests(shared_datadir):
        ds = stitch_dataset(
            counts=sim_test.pds.counts.to_pandas(),
            peptide_table=sim_test.pds.peptide_table.to_pandas(),
            sample_table=sim_test.pds.sample_table.to_pandas(),
        )
        assert sim_test.pds.equals(ds)
