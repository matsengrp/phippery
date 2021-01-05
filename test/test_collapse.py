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
from phippery.collapse import collapse_sample_groups


def test_collapse_sample_groups():

    fastq_filename = [f"sample_{i}.fastq" for i in range(12)]
    reference = [f"refa" for _ in range(12)]
    seq_dir = [f"expa" for _ in range(12)]
    bio_reps = [i for i in range(3) for _ in range(4)]
    tech_reps = [i for i in range(6) for _ in range(2)]

    columns = [
        "sample_id",
        "tech_rep_id",
        "bio_rep_id",
        "fastq_filename",
        "reference",
        "seq_dir",
    ]
    df = pd.DataFrame(
        zip(
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
    counts = np.ones([10, 12])
    ds = generate_sim_ds(counts=counts, sample_metadata=sample_metadata)
    # print(ds.counts)

    mean_tech_rep_ds = collapse_sample_groups(ds, group="tech_rep_id")
    # print(mean_tech_rep_ds.counts)

    expected_counts_solution = np.ones([10, 6])
    assert np.allclose(mean_tech_rep_ds.counts.values, expected_counts_solution)
    assert mean_tech_rep_ds.attrs["sample_coord_dim"] == "tech_rep_id"
