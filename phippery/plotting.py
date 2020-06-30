"""
@File: plotting.py

@Author: Jared Galloway

This should include some helpful function which
will allow us to plot interesting things given
a PhipData xarray object
"""


import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import collections as mc
from matplotlib import cm
from phippery.utils import get_std_dev_non_control, get_exp_prot_subset
import numpy as np
import xarray as xr
import pandas as pd
import itertools
import scipy.stats as st
import copy


def plot_tech_rep_corr(
    ds,
    color_by="sample_type",
    saveas=None,
    show=True,
    cmap="viridis",
    figsize=(15, 20),
    title_add="",
):
    """
    plot the technical replicate correlation for all samples in the xarray dataset
    provided. It is required that the "tech_rep_correlation" position exists
    in the coordinates. This function only plots samples that have 2 or more technical replicates

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: color_by <String> - The sample metadata ds you would like to color the bars by.
        This means there will be a unique color (according to cmap).

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.

    :param: cmap <String> - A string representing a matplotlib colormap.
    """

    cmap = plt.get_cmap(cmap)
    dss = {
        x for x in set(ds.sample_table.sel(sample_metadata=color_by).values) if x == x
    }
    colors = cmap(np.linspace(0, 1, len(dss)))
    color_dict = {ds: col for ds, col in zip(list(dss), colors)}
    black = [0.0, 0.0, 0.0, 1.0]
    color_dict["nan"] = black
    fig, ax = plt.subplots(1, figsize=figsize)
    sample_meta = ds.sample_table.to_pandas()
    two_replicate_samples = sample_meta[sample_meta["num_tech_reps"] == 2]
    control_st = list(two_replicate_samples.loc[:, f"{color_by}"].astype(str))
    sam_id = list(two_replicate_samples.index)
    corr = two_replicate_samples.tech_rep_correlation.values
    sort_a_by_b = lambda a, b: [  # noqa
        x for _, x in sorted(zip(b, a), key=lambda pair: pair[0])
    ]
    sorted_con_st = sort_a_by_b(control_st, corr)
    sorted_sam_id = [str(ID) for ID in sort_a_by_b(sam_id, corr)]
    sorted_corr = sorted(corr)
    sorted_colors = [color_dict[g] for g in sorted_con_st]
    ax.bar(sorted_sam_id, sorted_corr, color=sorted_colors)
    ax.set_title(f"{title_add}")
    ax.set_ylabel("Pearson Correlation")
    ax.set_xlabel("Sample ID")
    ax.yaxis.grid()
    markers = [
        plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
        for color in color_dict.values()
    ]
    ax.legend(
        markers, color_dict.keys(), bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0,
    )
    plt.tight_layout()
    if saveas:
        fig.savefig(f"{saveas}")  # noqa
    if show:
        plt.show()  # noqa


def plot_tech_rep_corr_split_by(
    ds,
    split_by="experiment",
    color_by="sample_type",
    saveas=None,
    show=True,
    cmap="viridis",
    figsize=(15, 20),
):
    """
    plot the technical replicate correlation for all samples in the xarray dataset
    provided. It is required that the "tech_rep_correlation" position exists
    in the coordinates. This function only plots samples that have 2 or more technical replicates

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: split_by <String> - The sample metadata group you would like to split subplots by.
        This means that there will be exactly one subplot for each unique factor.

    :param: color_by <String> - The sample metadata group you would like to color the bars by.
        This means there will be a unique color (according to cmap).

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.

    :param: cmap <String> - A string representing a matplotlib colormap.
    """

    cmap = plt.get_cmap(cmap)
    groups = {
        x for x in set(ds.sample_table.sel(sample_metadata=color_by).values) if x == x
    }
    colors = cmap(np.linspace(0, 1, len(groups)))
    color_dict = {group: col for group, col in zip(list(groups), colors)}
    black = [0.0, 0.0, 0.0, 1.0]
    color_dict["nan"] = black
    num_groups = len(set(ds.sample_table.loc[:, split_by].values))
    fig, ax = plt.subplots(num_groups, figsize=figsize)
    for idx, (experiment, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, split_by])
    ):
        sample_meta = group.sample_table.to_pandas()
        two_replicate_samples = sample_meta[sample_meta["num_tech_reps"] == 2]
        control_st = list(two_replicate_samples.loc[:, f"{color_by}"].astype(str))
        sam_id = list(two_replicate_samples.index)
        corr = two_replicate_samples.tech_rep_correlation.values
        sort_a_by_b = lambda a, b: [  # noqa
            x for _, x in sorted(zip(b, a), key=lambda pair: pair[0])
        ]
        sorted_con_st = sort_a_by_b(control_st, corr)
        sorted_sam_id = [str(ID) for ID in sort_a_by_b(sam_id, corr)]
        sorted_corr = sorted(corr)
        sorted_colors = [color_dict[g] for g in sorted_con_st]
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
    plt.tight_layout()
    if saveas:
        fig.savefig(f"{saveas}")  # noqa
    if show:
        plt.show()  # noqa


def plot_counts(
    ds,
    split_by="experiment",
    saveas=None,
    show=False,
    control_status_column="control_status",
    sample_label="empirical",
    log=True,
    figsize=(20, 20),
):
    """
    A Generalized function for plotting the counts table. This function
    requires that you have a sample metadata coordinate which deciphers
    between the controls and samples, we only plot samples here.

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: split_by <String> - The sample metadata column you would like to split by.
        This must exist in the sample_metadata coordinate of the dataset provided.

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.

    :param: control_status_column <String> - The name of the column to find control
        status of a sample.

    :param: sample_label <String> - The name of the factor in the control_status_column
        that represents the biological samples.

    :param: log <Bool> - Whether or not to log transform the matrix before plotting.

    :param: figsize <Tuple> - A Tuple specifying integer width and height of the the plot.
    """

    if control_status_column not in ds.sample_metadata.values:
        raise ValueError(
            f"{control_status_column} is not in sample_metadata coordinate"
        )
    num_groups = len(set(ds.sample_table.loc[:, split_by].values))
    fig, ax = plt.subplots(1, num_groups, figsize=figsize, sharey=True)
    for idx, (experiment, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, split_by])
    ):
        sam_meta_df = group.sample_table.to_pandas()
        emp_samples = sam_meta_df[
            sam_meta_df[control_status_column] == sample_label
        ].index
        empirical_counts = group.counts.loc[:, emp_samples]
        mat = np.log(empirical_counts + 0.1) if log else empirical_counts
        ax[idx].imshow(mat, aspect="auto")
        ax[idx].set_title(f"{experiment}")
        ax[idx].set_xlabel("Samples")
    ax[0].set_ylabel("Peptide ID")
    if saveas:
        fig.savefig(f"{saveas}")  # noqa
    if show:
        plt.show()  # noqa


def plot_enrichments(
    ds,
    exp_name,
    locus_name,
    split_by="reference",
    experiment_column="experiment",
    locus_name_column="Virus",
    locus_column="Prot_Start",
    saveas=None,
    show=True,
):
    """
    plot the counts matrix enrichment for a certain set of samples,
    split by grouping in the sample metadata table.

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: exp_name <String> - String representing the factor in
        'experiment_column' to subset the `ds` by before plotting.

    :param: split_by <String> - The sample metadata column you would like to split by.
        This must exist in the sample_metadata coordinate of the dataset provided.

    :param: experiment_column <String> the name of the column to find
        unique experiments in.

    :param: locus_name_column <String> the name of the locus column.

    :param: locus_column <String> the name of the columns specifying
        some information about the genomic location. The enrichment
        will be plotted as a function of this factor.

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.
    """

    if experiment_column not in ds.sample_metadata.values:
        raise ValueError(f"{experiment_column} is not in sample_metadata coordinate")

    if locus_name_column not in ds.peptide_metadata.values:
        raise ValueError(f"{locus_name_column} is not in peptide_metadata coordinate")

    if locus_column not in ds.peptide_metadata.values:
        raise ValueError(f"{locus_column} is not in peptide_metadata coordinate")

    exp_group = {
        exp: group
        for exp, group in ds.groupby(ds.sample_table.loc[:, experiment_column])
    }[exp_name]
    exp_strain_ds = {
        virus: group
        for virus, group in exp_group.groupby(
            exp_group.peptide_table.loc[:, locus_name_column]
        )
    }[locus_name]
    num_groups = len(
        [x for x in set(exp_strain_ds.sample_table.loc[:, split_by].values) if x == x]
    )
    fig, ax = plt.subplots(num_groups, figsize=(20, 20), sharey=True)
    ax = [ax] if num_groups == 1 else ax
    prot_start = exp_strain_ds.peptide_table.loc[:, locus_column]
    for idx, (label, group) in enumerate(
        exp_strain_ds.groupby(exp_strain_ds.sample_table.loc[:, split_by])
    ):
        ax[idx].set_title(
            f" experiment: {exp_name}\n protein: {locus_name}\n group: {label}"
        )
        for sample in group.sample_id:
            sample_enrichment = group.counts.loc[:, sample]
            ax[idx].plot(
                prot_start, sample_enrichment, linewidth=4, label=sample.values
            )
        ax[idx].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0)
    plt.subplots_adjust(hspace=0.4)
    if saveas:
        fig.savefig(f"{saveas}")  # noqa
    if show:
        plt.show()  # noqa


def biological_rep_correlation(
    ds,
    exp_name_list,
    color_by="reference",
    saveas=None,
    show=True,
    experiment_column="experiment",
    bio_id_column="sample_ID",
    cmap="viridis",
    title_add="",
    figsize=(10, 15),
):
    """
    plot the biological replicates for a set of experiments listed in
    'exp_name_list'. Each biological replicate label will also contain
    information about the samples and their respective technical replicate
    correlation.

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: experiment_name_list <List> - a list of strings representing the
        factors of 'experiment_column' to be included in the plot

    :param: color_by <String> - The sample_metadata coordinate grouping to color
        the bars by in the plot - a unique color pulled from `cmap` for each unique
        factor in the sample metadata coordinate.

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.

    :param: experiment_column <String> the name of the column
        to find unique experiments in.

    :param: bio_id_column <String> the column representing the biological
        id from which the sample itself derived. This is how we will group
        the biological replicates and calculate correlation.

    :param: cmap <String> - A string representing a matplotlib colormap.

    :param: title_add <String> - Anything extra you would like to add to the title
        of the plot.
    """

    if experiment_column not in ds.sample_metadata.values:
        raise ValueError("{experiment_column} does not exist in sample metadata")

    if bio_id_column not in ds.sample_metadata.values:
        raise ValueError("{bio_id_column} does not exist in sample metadata")

    # subset only the experiments we want.
    exp_group = {
        exp: group
        for exp, group in ds.groupby(ds.sample_table.loc[:, experiment_column])
        if exp in exp_name_list
    }
    exp_group = xr.combine_by_coords(
        list(exp_group.values()), data_vars="different", coords="different"
    )

    labels, pw_cc, col_group = [], [], []
    for emp_id, group in exp_group.groupby(
        exp_group.sample_table.loc[:, bio_id_column]
    ):
        col_group.append(group.sample_table.loc[:, color_by].values[0])
        bio_replicates = list(group.sample_id.values)
        respective_tech_rep_cc = group.sample_table.loc[
            bio_replicates, "tech_rep_correlation"
        ].values
        emp_sam_label = f"BIO ID: {emp_id} - REPLICATES: " + ", ".join(
            [
                f"{sam_id} : {round(tech_rep_cor, 2)}"
                for sam_id, tech_rep_cor in zip(bio_replicates, respective_tech_rep_cc)
            ]
        )
        labels.append(emp_sam_label)
        correlations = []
        for sample_ids in itertools.combinations(group.sample_id.values, 2):
            sample_0_enrichment = group.counts.loc[:, sample_ids[0]]
            sample_1_enrichment = group.counts.loc[:, sample_ids[1]]
            correlation = (
                st.pearsonr(sample_0_enrichment, sample_1_enrichment)[0]
                if np.any(sample_0_enrichment != sample_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(sum(correlations) / len(correlations))

    cmap = plt.get_cmap(cmap)
    groups = {x for x in set(col_group) if x == x}
    colors = cmap(np.linspace(0, 1, len(groups)))
    color_dict = {group: col for group, col in zip(list(groups), colors)}
    black = [0.0, 0.0, 0.0, 1.0]
    color_dict[np.nan] = black

    fig, ax = plt.subplots(1, figsize=(10, 10))
    sort_a_by_b = lambda a, b: [x for _, x in sorted(zip(b, a))]  # noqa
    sorted_labels = sort_a_by_b(labels, pw_cc)
    sorted_col_group = sort_a_by_b(col_group, pw_cc)
    sorted_corr = sorted(pw_cc)
    sorted_colors = []
    # is there a cleaner way to deal with NaN's?
    for con_st in sorted_col_group:
        if con_st != con_st:
            sorted_colors.append(black)
        else:
            sorted_colors.append(color_dict[con_st])
    ax.barh(np.arange(len(sorted_labels)), sorted_corr, color=sorted_colors)
    ax.set_yticks(np.arange(len(sorted_labels)))
    ax.set_yticklabels(sorted_labels)
    exps_str = " | ".join(exp_name_list)
    ax.set_title(f"Biological Replicate Correlation \n Exp: {exps_str}\n{title_add}")
    ax.set_xlabel("Pearson Correlation")
    ax.set_ylabel("Empirical Sample ID")
    ax.xaxis.grid()
    markers = [
        plt.Line2D([0, 0], [0, 0], color=color, marker="o", linestyle="")
        for color in color_dict.values()
    ]
    ax.legend(
        markers, color_dict.keys(), bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0
    )
    plt.tight_layout()
    if saveas:
        fig.savefig(f"{saveas}")  # noqa
    if show:
        plt.show()  # noqa


def plot_temporal_enrichments(
    ds,
    temporal_name_column="days_from_symptom_onset",
    saveas=None,
    show=True,
    figsize=(10, 15),
    title_add="",
):
    """
    A general function for plotting a heatmap showing the
    mean enrichment for a set of samples grouped by a the temporal
    column provided.

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: temporal_name_column <String> - the name of the column have numerical
        temporal data.

    :param: saveas <String> - The path you would like to save the plot to.

    :param: show <Bool> - Whether or not to use matplotlib show() function.

    :param: figsize <Tuple> - figure size.

    :param: title_add <String> - Any string you would like to add to the title.
    """

    if temporal_name_column not in ds.sample_metadata.values:
        raise ValueError("{temporal_name_column} does not exist in sample metadata")

    sd = round(get_std_dev_non_control(ds), 2)
    all_unique_times = [
        x for x in set(ds.sample_table.loc[:, temporal_name_column].values) if x == x
    ]
    temporal_enrichments = np.zeros((len(all_unique_times), len(ds.peptide_id.values)))
    times, num_samples_per_group = [], []
    for time_index, (time, group) in enumerate(
        ds.groupby(ds.sample_table.loc[:, temporal_name_column])
    ):
        times.append(time)
        num_samples_per_group.append(len(group.sample_id.values))
        temporal_enrichments[time_index, :] += group.counts.values.sum(axis=1)

    nte = temporal_enrichments / np.array(num_samples_per_group)[:, None]
    df = pd.DataFrame(nte, columns=ds.peptide_id.values, index=all_unique_times)
    df.sort_index(inplace=True)
    fig, ax = plt.subplots(1, figsize=figsize, sharey=True)
    sns.heatmap(df, ax=ax)
    ax.set_title(
        f"Mean Enrichment\n non control sample $\sigma$ = {sd}\n{title_add}"  # noqa
    )
    ax.set_ylabel("Days from symptom onset")
    ax.set_xlabel("Peptide id")
    if show:
        plt.show()
    if saveas:
        fig.savefig(saveas)
