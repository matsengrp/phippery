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
from phippery.collapse import pairwise_correlation_by_sample_group
from phippery.normalize import cpm


def test_collapse_sample_groups():

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
    counts = np.random.randint(0, 10, [10, 12])
    ds = generate_sim_ds(counts=counts, sample_metadata=sample_metadata)
    cpm(ds, per_sample=True, inplace=True)

    mean_tech_rep_ds = collapse_sample_groups(
        ds, group="tech_rep_id", compute_pw_cc=True
    )

    # FIXME This should do it -> other than the datatype issue of counts
    # mean_tech_rep_ds = collapse_sample_groups(
    #    ds, group="u_id", compute_pw_cc=True
    # )

    assert "tech_rep_id_cpm_pw_cc" in mean_tech_rep_ds.sample_metadata.values
    assert "tech_rep_id_cpm_n_reps" in mean_tech_rep_ds.sample_metadata.values
    assert "tech_rep_id_counts_pw_cc" in mean_tech_rep_ds.sample_metadata.values
    assert "tech_rep_id_counts_n_reps" in mean_tech_rep_ds.sample_metadata.values

    cpm_ex = ds["cpm"].values
    comb = [(i, i + 1) for i in range(0, 12, 2)]

    for m in ["counts", "cpm"]:

        enr = counts if m == "counts" else cpm_ex

        mean_sol = np.zeros([10, 6])
        for t, (i, j) in enumerate(comb):
            mean_sol[:, t] = np.mean(enr[:, [i, j]], axis=1)
        local_mean_sol = mean_tech_rep_ds[m].values

        # compare our solution to a more interpretable solution
        assert np.allclose(mean_sol, local_mean_sol)

        enr_corr = pd.DataFrame(enr).corr()
        cc_sol = np.array([round(enr_corr.iloc[i, j], 5) for (i, j) in comb])
        local_cc_sol = mean_tech_rep_ds.sample_table.loc[
            :, f"tech_rep_id_{m}_pw_cc"
        ].values.astype(np.float64)

        # compare our solution to pandas corr computed manually
        assert np.allclose(cc_sol, local_cc_sol)


def test_pairwise_correlation_by_sample_group():
    # TODO
    pass
