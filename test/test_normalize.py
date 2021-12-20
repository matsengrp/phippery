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
from phippery.normalize import _comp_std_enr
from phippery.normalize import _wt_window_scalar
from phippery.normalize import _comp_enr
from phippery.normalize import _comp_diff_sel
from phippery.normalize import _comp_size_factors
from phippery.normalize import size_factors
from phippery.normalize import _comp_cpm
from phippery.normalize import counts_per_million
from phippery.normalize import _comp_cpm_per_sample
from phippery.normalize import _comp_rank
from phippery.normalize import rank_data
from phippery.normalize import _comp_rank_per_sample
from phippery.normalize import differential_selection_wt_mut
from phippery.normalize import differential_selection_sample_groups
from phippery.normalize import replicate_oligo_counts
from sim_test_generator import generate_sim_ds
from sim_test_generator import make_hardcoded_ds


def test_comp_enr():
    """
    test the enrichment of all ones
    """

    for i in range(3, 7):
        enr = _comp_enr(
            counts=pd.DataFrame(np.ones([i, 2])).astype(float),
            lib_counts=pd.DataFrame(np.ones([i, 1])).astype(float),
        )
        assert np.all(enr == np.ones([i, 2]))


def test_comp_std_enr():
    """
    test the standard enrichment of all ones which should result in zeros in all places
    """

    for i in range(3, 7):
        std_enr = _comp_std_enr(
            counts=pd.DataFrame(np.ones([i, 2])).astype(float),
            lib_counts=pd.DataFrame(np.ones([i, 1])).astype(float),
            mock_ip_counts=pd.DataFrame(np.ones([i, 1])).astype(float),
        )
        assert np.all(std_enr == np.zeros([i, 2]))


def test_comp_diff_sel():

    assert np.all(_comp_diff_sel(1, [2, 2, 2], scalar=1) == [1.0, 1.0, 1.0])
    assert np.all(_comp_diff_sel(2, [1, 1, 1], scalar=1) == [-1.0, -1.0, -1.0])
    assert np.all(_comp_diff_sel(2, [1, 1, 1], scalar=2) == [-2.0, -2.0, -2.0])
    assert np.all(_comp_diff_sel(1, [1, 1, 1], scalar=1) == [0, 0, 0])
    assert np.all(_comp_diff_sel(1, []) == [])


def test_zero_div_error():

    with pytest.raises(ZeroDivisionError):
        _comp_diff_sel(1, [0])


def test_wt_scalar():

    assert _wt_window_scalar([2, 0, 1, 2, 0], 0, 2) == 1.0
    assert _wt_window_scalar([2, 0, 1, 2, 0], 1, 2) == 1.25
    assert _wt_window_scalar([2, 0, 1, 2, 0], 2, 2) == 1.0
    assert _wt_window_scalar([2, 0, 1, 2, 0], 3, 2) == 0.75
    assert _wt_window_scalar([2, 0, 1, 2, 0], 4, 2) == 1.0


def test_diff_sel_wt_mut():

    prot = [f"prot_{l}" for l in ["a", "b"] for _ in range(100)]
    loc = [i for i in range(10) for _ in range(10)] * 2
    is_wt = ([True] + ([False] * 9)) * 20
    oligo = ["ATCG"] * 200

    pep_df = pd.DataFrame(
        zip(range(200), oligo, prot, loc, is_wt),
        columns=["peptide_id", "oligo", "prot", "loc", "is_wt"],
    )

    counts = np.full([200, 10], 2)
    for row in range(len(pep_df)):
        if pep_df.loc[row, "is_wt"]:
            counts[row, :] = 1

    ds = generate_sim_ds(counts=counts, sample_metadata=None, peptide_metadata=pep_df)

    differential_selection_wt_mut(
        ds,
        data_table="counts",
        scaled_by_wt=True,
        smoothing_flank_size=5,
        groupby=["prot"],
        loc_column="loc",
        is_wt_column="is_wt",
        inplace=True,
        new_table_name="ds",
    )

    counts_copy = copy.deepcopy(counts) - 1
    assert np.allclose(counts_copy, ds["ds"].values)


def test_svd_aa_loc():

    pass


def test_differential_selection_sample_groups():

    num_samples = 4
    fastq_filename = [f"sample_{i}.fastq" for i in range(num_samples)]
    reference = [f"refa" for _ in range(num_samples)]
    seq_dir = [f"expa" for _ in range(num_samples)]
    library_batch = [f"batch_{l}" for l in ["a", "b"] for _ in range(num_samples // 2)]

    columns = ["sample_id", "fastq_filename", "reference", "seq_dir", "library_batch"]
    df = pd.DataFrame(
        zip(range(num_samples), fastq_filename, reference, seq_dir, library_batch),
        columns=columns,
    )
    sample_metadata = df.set_index("sample_id")

    counts = np.full([5, 4], 2)
    counts[:, [0, 1]] = 1

    ds = generate_sim_ds(counts=counts, sample_metadata=sample_metadata)

    differential_selection_sample_groups(
        ds,
        sample_feature="library_batch",
        is_equal_to="batch_a",
        data_table="counts",
        aggregate_function=np.mean,
        inplace=True,
        new_table_name="ds",
    )
    counts_copy = copy.deepcopy(counts) - 1
    assert np.allclose(counts_copy, ds["ds"].values)


def test_size_factors():
    # fix the masking error then head here.
    """
    A single masked value should not effect the values being normalized if all
    of the values are equal to 1
    """
    arr = np.ones([6, 4])
    arr[0, 0] = 0
    sf = _comp_size_factors(arr)
    assert np.allclose(sf, arr)


def test_size_factors_ones():
    """
    we expect all ones to go to infinity.
    """

    arr = np.ones([10, 16])
    sf = _comp_size_factors(arr)
    assert np.all(sf == np.full([10, 16], np.inf))


def test_size_factors_regular():

    ds = generate_sim_ds()
    size_factors(ds)
    assert "size_factors" in ds.data_vars


def test_cpm():

    for wh in range(2, 10):
        arr = np.ones([wh, wh])
        cpm = _comp_cpm(arr)
        expectation = (1 / (wh * wh)) * 1e6
        assert np.allclose(cpm, expectation)


def test_cpm_per_sample():

    for wh in range(2, 10):
        arr = np.ones([wh, wh])
        cpm = _comp_cpm_per_sample(arr)
        expectation = (1 / wh) * 1e6
        assert np.allclose(cpm, expectation)


# TODO add this after implimenting simple sim generator
def test_cpm_ds():

    ds = generate_sim_ds()
    counts_per_million(ds)
    assert "cpm" in ds.data_vars


def test_rank_data():

    arr = np.array([[0, 1], [2, 3]])
    rank = _comp_rank(arr)
    assert np.all(arr == rank)


def test_rank_data_per_sample():

    arr = np.repeat(range(10), 5).reshape(10, 5)
    rank = _comp_rank_per_sample(arr)
    assert np.all(arr == rank)


def test_rank_ds():

    ds = generate_sim_ds()
    rank_data(ds)
    assert "rank" in ds.data_vars


def test_replicate_oligo_counts(shared_datadir):
    """
    assert that replicates peptides in the library have equal counts
    after using the utils.replicate_oligo_counts function.
    """
    ds = make_hardcoded_ds()
    peptide_sums = ds.counts.sum(axis=0).values
    replicate_oligo_counts(ds)
    counts = ds.counts.to_pandas()
    for i, (sample, pep_enr) in enumerate(counts.iteritems()):
        assert np.all(pep_enr == peptide_sums[i])
