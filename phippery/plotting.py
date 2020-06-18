"""
@File: plotting.py

@Author: Jared Galloway

This should include some helpful function which
will allow us to plot interesting things given
a PhipData object
"""

import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import cm
from phippery.PhipData import *
import numpy as np


def plot_tech_rep_corr_xarray(
    ds, split_by="experiment", color_by="sample_type", saveas=None, show=True
):
    cmap = plt.get_cmap("viridis")
    groups = {
        x for x in set(ds.sample_table.sel(sample_metadata=color_by).values) if x == x
    }
    colors = cmap(np.linspace(0, 1, len(groups)))
    color_dict = {group: col for group, col in zip(list(groups), colors)}
    black = [0.0, 0.0, 0.0, 1.0]
    color_dict[np.nan] = black
    num_groups = len(set(ds.sample_table.loc[:, split_by].values))
    fig, ax = plt.subplots(num_groups, figsize=(20, 25))
    for idx, (experiment, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, split_by])
    ):
        sample_meta = group.sample_table.to_pandas()
        two_replicate_samples = sample_meta[sample_meta["num_tech_reps"] == 2]
        control_st = two_replicate_samples[f"{color_by}"]
        sam_id = two_replicate_samples.index
        corr = two_replicate_samples.tech_rep_correlation
        sort_a_by_b = lambda a, b: [x for _, x in sorted(zip(b, a))]  # noqa
        sorted_con_st = sort_a_by_b(control_st, corr)
        sorted_sam_id = [str(ID) for ID in sort_a_by_b(sam_id, corr)]
        sorted_corr = sorted(corr)
        sorted_colors = []
        # is there a cleaner way to deal with NaN's?
        for con_st in sorted_con_st:
            if con_st != con_st:
                sorted_colors.append(black)
            else:
                sorted_colors.append(color_dict[con_st])
        ax[idx].bar(sorted_sam_id, sorted_corr, color=sorted_colors)
        ax[idx].set_title(f"{experiment}")
        ax[idx].set_ylabel("Pearson Correlation")
        ax[idx].set_xlabel("Sample ID")
        ax[idx].yaxis.grid()
        markers = [
            plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
            for color in color_dict.values()
        ]
        ax[idx].legend(
            markers,
            color_dict.keys(),
            bbox_to_anchor=(1.01, 1),
            loc=2,
            borderaxespad=0.0,
        )
    plt.subplots_adjust(hspace=0.5)
    if saveas:
        fig.savefig(f"{saveas}")
    if show:
        plt.show()


def plot_counts(ds, split_by="experiment", saveas=None, show=True):

    assert "control_status" in ds.sample_metadata.values
    num_groups = len(set(ds.sample_table.loc[:, split_by].values))
    fig, ax = plt.subplots(1, num_groups, figsize=(20, 20), sharey=True)
    for idx, (experiment, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, split_by])
    ):
        sam_meta_df = group.sample_table.to_pandas()
        emp_samples = sam_meta_df[sam_meta_df["control_status"] == "empirical"].index
        empirical_counts = group.counts.loc[:, emp_samples]
        empirical_counts += 0.1  # so we can do log
        ax[idx].imshow(np.log(empirical_counts), aspect="auto")
        ax[idx].set_title(f"{experiment}")
        ax[idx].set_xlabel("Samples")
    ax[0].set_ylabel("Peptide ID")
    if saveas:
        fig.savefig(f"{saveas}")
    if show:
        plt.show()


def compute_standardized_enrichment(ds):
    """
    return a new xarray dataset same as the input
    except with the counts converted to standard enrichment.

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5
    """

    assert "experiment" in ds.sample_metadata.values
    ret = ds.copy()

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


def plot_enrichments(
    ds,
    exp_name,
    protein_name,
    sample_split_column="reference",
    pseudo_count_bias=5,
    saveas=None,
    show=True,
):
    assert (
        "experiment" in ds.sample_metadata.values
        and "Virus" in ds.peptide_metadata.values
    )
    # TODO now assert the actual experiment and virus exist!
    exp_group = {
        exp: group for exp, group in ds.groupby(ds.sample_table.loc[:, "experiment"])
    }[exp_name]
    exp_strain_ds = {
        virus: group
        for virus, group in exp_group.groupby(exp_group.peptide_table.loc[:, "Virus"])
    }[protein_name]
    num_groups = len(
        [
            x
            for x in set(exp_strain_ds.sample_table.loc[:, sample_split_column].values)
            if x == x
        ]
    )
    fig, ax = plt.subplots(num_groups, figsize=(20, 20), sharey=True)
    ax = [ax] if num_groups == 1 else ax
    prot_start = exp_strain_ds.peptide_table.loc[:, "Prot_Start"]
    for idx, (label, group) in enumerate(
        exp_strain_ds.groupby(exp_strain_ds.sample_table.loc[:, sample_split_column])
    ):
        ax[idx].set_title(
            f" experiment: {exp_name}\n protein: {protein_name}\n group: {label}"
        )
        for sample in group.sample_id:
            sample_enrichment = group.counts.loc[:, sample]
            ax[idx].plot(
                prot_start, sample_enrichment, linewidth=4, label=sample.values
            )
        ax[idx].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
    plt.subplots_adjust(hspace=0.4)
    if saveas:
        fig.savefig(f"{saveas}")
    plt.show()
