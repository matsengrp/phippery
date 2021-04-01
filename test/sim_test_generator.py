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

from phippery.phipdata import stitch_dataset

simulation_tests = ["simulate_small_ones_sep_reps"]


# TODO add more helpful attributes
class SimulationTest(object):

    """
    A class that can define the aspects
    of each simulation test.
    """

    def __init__(self, path):

        xr_pd = os.path.join(path, "phip_data/pds.phip")
        self.pds = pickle.load(open(xr_pd, "rb"))
        sol = os.path.join(path, "solution.np")
        self.solution = pickle.load(open(sol, "rb"))
        self.counts = [f for f in glob.glob(os.path.join(path, "counts/*/*"))]
        self.sam_meta = os.path.join(path, "sample_metadata.csv")
        self.pep_meta = os.path.join(path, "peptide_metadata.csv")
        self.name = str(os.path.dirname(path))


def iter_sim_tests(shared_datadir):
    """
    a simple generator for the existing
    simulation tests.
    """

    for test in simulation_tests:
        assert os.path.exists((shared_datadir / test))
        yield SimulationTest((shared_datadir / test))


def generate_sim_ds(
    counts=np.ones([5, 4]), sample_metadata=None, peptide_metadata=None
):
    """
    Generate a simple simulated dataset with the bare minimum
    unless provided with more complex tables.
    """

    if sample_metadata is None:
        num_samples = counts.shape[1]
        fastq_filename = [f"sample_{i}.fastq" for i in range(num_samples)]
        reference = [f"refa" for _ in range(num_samples)]
        seq_dir = [f"expa" for _ in range(num_samples)]

        columns = ["sample_id", "fastq_filename", "reference", "seq_dir"]
        df = pd.DataFrame(
            zip(range(num_samples), fastq_filename, reference, seq_dir), columns=columns
        )
        sample_metadata = df.set_index("sample_id")

    if peptide_metadata is None:
        num_peptides = counts.shape[0]
        columns = ["peptide_id", "Oligo"]
        oligos = ["ATCG" for _ in range(num_peptides)]
        df = pd.DataFrame(zip(range(num_peptides), oligos), columns=columns)
        peptide_metadata = df.set_index("peptide_id")

    counts_df = pd.DataFrame(counts)
    return stitch_dataset(counts_df, peptide_metadata, sample_metadata)
