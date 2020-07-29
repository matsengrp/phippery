"""
@File: sim_test_generator.py

@Author: Jared Galloway

A function which can generate simulation
class tests.
"""

# built-in dependency
import os
import sys

# dependency
import pytest
import numpy as np
import pandas as pd
import xarray as xr
import pickle
import glob


simulation_tests = ["simulate_small_ones", "simulate_medium_1_mismatch"]


# TODO add more helpful attributes
class SimulationTest(object):

    """
    A class that can define the aspects
    of each simulation test.
    """

    def __init__(self, path):

        xr_pd = os.path.join(path, "data/phip_data/pds.phip")
        self.pds = pickle.load(open(xr_pd, "rb"))
        sol = os.path.join(path, "solution.np")
        self.solution = pickle.load(open(sol, "rb"))

    # TODO
    def get_counts_pds_stub(self):
        pass


def iter_sim_tests(shared_datadir):
    """
    a simple generator for the existing
    simulation tests.
    """

    for test in simulation_tests:
        assert os.path.exists((shared_datadir / test))
        yield SimulationTest((shared_datadir / test))
