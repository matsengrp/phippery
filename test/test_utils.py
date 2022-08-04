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
from sim_test_generator import generate_sim_ds

# functions to test
from phippery.utils import to_tall
from phippery.utils import id_query
from phippery.utils import ds_query
from phippery.utils import get_annotation_table

from phippery.utils import to_wide_csv
from phippery.utils import dataset_from_csv


def test_tall_ds_shape():

    ds = generate_sim_ds()
    tall = to_tall(ds)
    num_counts = ds.counts.shape[0] * ds.counts.shape[1]
    assert num_counts == len(tall)


def test_query():

    ds = make_hardcoded_ds()
    wt_idx = id_query(ds, "is_wt == True", dim="peptide")
    assert np.all(wt_idx == [0, 5])
    ds_slice = ds_query(ds, "participant_id == 1")
    assert np.all(ds_slice.sample_id.values == [4, 5, 6, 7])

#def test_create_file(tmp_path)
#    d = tmp_path / "sub"
#    d.mkdir()
#    p = d / "hello.txt"
#    p.write_text(CONTENT)


def test_csv_functions(tmp_path):
    """assert that the to_csv is the inverse of dataset_from_csv"""

    d = tmp_path / "sub"
    d.mkdir()
    ds_from_mem = make_hardcoded_ds()
    prefix = "test"
    to_wide_csv(ds_from_mem, d / prefix)
    ds_from_disk = dataset_from_csv(
        peptide_table_filename = f"{d}/{prefix}_peptide_annotation_table.csv",
        sample_table_filename = f"{d}/{prefix}_sample_annotation_table.csv",
        counts_table_filename = f"{d}/{prefix}_counts.csv"
    )
    assert ds_from_disk.equals(ds_from_mem)
