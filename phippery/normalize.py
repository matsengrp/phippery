"""
@File: normalize.py

@Author: Jared Galloway

this file will contain functions which take
an xarray phip dataset and return a copy
where the counts have been normalized or transformed.
"""

import numpy as np
import xarray as xr
import itertools
import copy


def compute_standardized_enrichment(ds):
    """
    return a new xarray dataset same as the input
    except with the counts converted to standard enrichment.

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        proveded to phippery.collect
    """

    assert "experiment" in ds.sample_metadata.values
    ret = copy.deepcopy(ds)

    # each enrichment is specific to experiment control
    for idx, (experiment, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, "experiment"])
    ):
        group_counts = group.counts.to_pandas()
        group_sample_meta = group.sample_table.to_pandas()

        # find controls and average all
        group_lib_control_indices = group_sample_meta[
            group_sample_meta["control_status"] == "library"
        ].index
        group_bead_control_indices = group_sample_meta[
            group_sample_meta["control_status"] == "beads_only"
        ].index
        assert (
            len(group_lib_control_indices) > 0 and len(group_bead_control_indices) > 0
        )
        group_lib_counts_mean = group_counts[group_lib_control_indices].mean(axis=1)
        group_bead_counts_mean = group_counts[group_bead_control_indices].mean(axis=1)
        group_lib_counts_mean_sum = sum(group_lib_counts_mean)

        # compute beads control std enrichment
        pseudo_sample = group_bead_counts_mean + max(
            1, sum(group_bead_counts_mean) / group_lib_counts_mean_sum
        )
        pseudo_lib_control = group_lib_counts_mean + max(
            1, group_lib_counts_mean_sum / sum(group_bead_counts_mean)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        pseudo_bead_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        for bead_id in group_bead_control_indices:
            ret.counts.loc[:, bead_id] = pseudo_bead_enrichment

        # compute all sample standardized enrichment
        for sample_id, sample in group_counts.iteritems():
            if sample_id in list(group_lib_control_indices) + list(
                group_bead_control_indices
            ):
                continue
            pseudo_sample = sample + max(1, sum(sample) / group_lib_counts_mean_sum)
            pseudo_lib_control = group_lib_counts_mean + max(
                1, group_lib_counts_mean_sum / sum(sample)
            )
            pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
            pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
            sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
            ret.counts.loc[:, sample_id] = sample_enrichment - pseudo_bead_enrichment

    return ret
