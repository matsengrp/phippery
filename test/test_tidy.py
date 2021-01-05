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

# functions I'll be testing here
from phippery.tidy import tidy_ds


def test_tidy_ds_shape():

    ds = generate_sim_ds()
    tidy = tidy_ds(ds)
    num_counts = ds.counts.shape[0] * ds.counts.shape[1]
    assert num_counts == len(tidy)


# TODO, write some more detailed tests
