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


def make_hardcoded_ds():

    """
    A hard coded annotation dataset for testing all sorts of
    functions

    COUNTS

    sample_id   0   1   2   3   4   5   6   7   8   9   10  11
    peptide_id
    0            3   6   8   9   6   8   7   9   3   6   1   2
    1            5   5   0   5   0   9   9   3   1   7   4   1
    2            1   4   6   7   3   6   9   2   3   0   8   6
    3            6   0   5   6   0   2   0   3   0   6   2   2
    4            5   7   9   4   7   7   3   5   5   6   4   6
    5            8   5   0   4   1   9   4   5   9   8   7   1
    6            0   5   3   7   6   4   0   2   5   1   8   2
    7            4   0   4   3   4   0   6   8   6   1   9   9
    8            5   6   5   9   5   7   4   0   6   9   8   0
    9            1   6   2   0   0   0   2   2   6   2   8   9

    SAMPLE ANNOTATIONS

    sample_metadata u_id tech_rep_id bio_rep_id   fastq_filename reference seq_dir
    sample_id
    0                  0           0          0   sample_0.fastq      refa    expa
    1                  1           0          0   sample_1.fastq      refa    expa
    2                  2           1          0   sample_2.fastq      refa    expa
    3                  3           1          0   sample_3.fastq      refa    expa
    4                  4           2          1   sample_4.fastq      refa    expa
    5                  5           2          1   sample_5.fastq      refa    expa
    6                  6           3          1   sample_6.fastq      refa    expa
    7                  7           3          1   sample_7.fastq      refa    expa
    8                  8           4          2   sample_8.fastq      refa    expa
    9                  9           4          2   sample_9.fastq      refa    expa
    10                10           5          2  sample_10.fastq      refa    expa
    11                11           5          2  sample_11.fastq      refa    expa

    PEPTIDE ANNOTATIONS

    peptide_metadata Oligo  is_wt
    peptide_id
    0                 ATCG  False
    1                 ATCG   True
    2                 ATCG   True
    3                 ATCG   True
    4                 ATCG   True
    5                 ATCG  False
    6                 ATCG   True
    7                 ATCG   True
    8                 ATCG   True
    9                 ATCG   True

    """

    # TODO test a non-collapse
    fastq_filename = [f"sample_{i}.fastq" for i in range(12)]
    reference = [f"refa" for _ in range(12)]
    seq_dir = [f"expa" for _ in range(12)]
    bio_reps = [i for i in range(3) for _ in range(4)]
    tech_reps = [i for i in range(6) for _ in range(2)]

    columns = [
        "sample_id",
        "u_id",
        "tech_rep_id",
        "bio_rep_id",
        "fastq_filename",
        "reference",
        "seq_dir",
    ]
    df = pd.DataFrame(
        zip(
            range(len(reference)),
            range(len(reference)),
            tech_reps,
            bio_reps,
            fastq_filename,
            reference,
            seq_dir,
        ),
        columns=columns,
    )
    sample_metadata = df.set_index("sample_id")

    columns = ["peptide_id", "Oligo", "is_wt"]
    oligos = ["ATCG" for _ in range(10)]
    is_wt = [True] * 10
    is_wt[0] = False
    is_wt[5] = False
    df = pd.DataFrame(zip(range(10), oligos, is_wt), columns=columns)
    peptide_metadata = df.set_index("peptide_id")

    np.random.seed(23)
    counts = np.random.randint(0, 10, [10, 12])

    ds = generate_sim_ds(
        counts=counts,
        sample_metadata=sample_metadata,
        peptide_metadata=peptide_metadata,
    )

    return ds
