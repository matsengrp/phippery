"""
@File: test_phip_data.py

@Author: Jared Galloway

A set of unit tests for the gathering and exporting of phip data.
"""

# built-in dependency
import os
import sys
import copy

# dependency
import pytest
import numpy as np
import pandas as pd
import xarray as xr
import glob

# local
from sim_test_generator import generate_sim_ds
from sim_test_generator import make_hardcoded_ds

# functions I'll be testing here
from phippery.utils import collapse_groups
from phippery.utils import _mean_pw_cc_by

from phippery.normalize import counts_per_million
from phippery.utils import stitch_dataset

def test_throw_mm_features():
    pass


def test_mean_pw_cc_peptide_runs():
    ds = make_hardcoded_ds()
    _mean_pw_cc_by(ds, by="is_wt", dim="peptide")
    pass


def test_mean_pw_cc_sample_runs():
    ds = make_hardcoded_ds()
    _mean_pw_cc_by(ds, by="participant_id")
    pass


def test_simgle_anno_collapse_peptides():

    ds = make_hardcoded_ds()
    dsb = collapse_groups(ds, by=["is_wt"], collapse_dim="peptide")

    counts = ds.counts.values
    new = np.zeros([2, 12])
    new[1, :] = counts[[0, 5], :].mean(axis=0)
    new[0, :] = counts[[1, 2, 3, 4, 6, 7, 8, 9], :].mean(axis=0)

    assert np.all(dsb.counts.values == new)


def test_single_anno_collapse_samples():

    ds = make_hardcoded_ds()
    dsb = collapse_groups(ds, by=["participant_id"], collapse_dim="sample")

    counts = ds.counts.values
    new = np.zeros([10, 3])
    new[:, 0] = counts[:, 0:4].mean(axis=1)
    new[:, 1] = counts[:, 4:8].mean(axis=1)
    new[:, 2] = counts[:, 8:12].mean(axis=1)

    assert np.all(dsb.counts.values == new)
