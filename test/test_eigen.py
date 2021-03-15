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

# functions I'll be testing here
from phippery.eigen import eigenassay_projections
from sim_test_generator import generate_sim_ds


def test_eigenassay_projections():

    """
    simple test of the projections code on a simulated ds
    """

    ds = generate_sim_ds()
    ep = eigenassay_projections(
        ds,
        compute_correlations=True,
        return_raw_decomposition=True,
        return_eigenassay_meta=True,
    )

    assert type(ep) == dict
    assert "sample_eigenassay_projections" in ep
