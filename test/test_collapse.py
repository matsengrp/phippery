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
from phippery.collapse import collapse_groups
from phippery.collapse import mean_pw_cc_by

# from phippery.collapse import pairwise_correlation_by_sample_group
from phippery.normalize import cpm

from phippery.phipdata import stitch_dataset

# TODO more edge cases
# More error handle testing.


# TODO
def test_throw_mm_features():
    pass


def test_mean_pw_cc_peptide_runs():
    ds = make_hardcoded_ds()
    mean_pw_cc_by(ds, by="is_wt", dim="peptide")
    pass


def test_mean_pw_cc_sample_runs():
    ds = make_hardcoded_ds()
    mean_pw_cc_by(ds, by="participant_id")
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
