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


def standardized_enrichment_across_exp(ds):
    """
    return a new xarray dataset same as the input
    except with the counts converted to standard enrichment.
    This method requires

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        proveded to phippery.collect
    """
    print("WARNING: DEPRECATED")

    # we are returning a completely new dataset.
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
        if len(group_lib_control_indices) == 0:
            raise ValueError(
                "Experiment {experiment} does not have an associated library input control."
            )
        if len(group_bead_control_indices) == 0:
            raise ValueError(
                "Experiment {experiment} does not have an associated beads only control."
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


def standardized_enrichment(ds, ds_lib_control_indices, ds_bead_control_indices):
    """
    return a new xarray dataset same as the input
    except with the counts converted to standard enrichment.

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: ds_lib_control_indices <list> - a list of integers specifying the
        sample id's of the library controls you would like to you normalize
        all other samples with. We take the average of all lib controls.

    :param: ds_bead_control_indices <list> - a list of integers specifying the
        sample id's of the bead controls you would like to you normalize
        all other samples with. We take the average of all lib controls.
    """

    # we are returning a completely new dataset.
    ret = copy.deepcopy(ds)

    ds_counts = ds.counts.to_pandas()

    # find controls and average all
    ds_lib_counts_mean = ds_counts[ds_lib_control_indices].mean(axis=1)
    ds_bead_counts_mean = ds_counts[ds_bead_control_indices].mean(axis=1)
    ds_lib_counts_mean_sum = sum(ds_lib_counts_mean)

    # compute beads control std enrichment
    pseudo_sample = ds_bead_counts_mean + max(
        1, sum(ds_bead_counts_mean) / ds_lib_counts_mean_sum
    )
    pseudo_lib_control = ds_lib_counts_mean + max(
        1, ds_lib_counts_mean_sum / sum(ds_bead_counts_mean)
    )
    pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
    pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
    pseudo_bead_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
    for bead_id in ds_bead_control_indices:
        ret.counts.loc[:, bead_id] = pseudo_bead_enrichment

    # compute all sample standardized enrichment
    for sample_id, sample in ds_counts.iteritems():
        if sample_id in list(ds_lib_control_indices) + list(ds_bead_control_indices):
            continue
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_control = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        ret.counts.loc[:, sample_id] = sample_enrichment - pseudo_bead_enrichment

    return ret


def enrichment(ds, ds_lib_control_indices):
    """
    return a new xarray dataset same as the input
    except with the counts converted to enrichment.

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :param: ds_lib_control_indices <list> - a list of integers specifying the
        sample id's of the library controls you would like to you normalize
        all other samples with. We take the average of all lib controls.
    """

    # we are returning a completely new dataset.
    ret = copy.deepcopy(ds)
    # ret = ds.copy(deep=True)

    ds_counts = ds.counts.to_pandas()

    # find controls and average all
    ds_lib_counts_mean = ds_counts[ds_lib_control_indices].mean(axis=1)
    ds_lib_counts_mean_sum = sum(ds_lib_counts_mean)

    # compute all sample standardized enrichment
    for sample_id, sample in ds_counts.iteritems():
        if sample_id in list(ds_lib_control_indices):
            continue
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_control = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        ret.counts.loc[:, sample_id] = sample_enrichment

    return ret


def differential_selection(
    ds,
    protein_name_column="Protein",
    wd_location_column="Loc",
    aa_sub_column="aa_sub",
    is_wt_column="is_wt",
):
    """
    Compute differential selection for all counts in a ds.

    Identifying the wild type at each location for each protein,
    This function computes log fold change between the wild type and
    mutant for each protein/loc combo.

    :param: ds <xarray.Dataset> - An xarray dataset obtained from three tables
        provided to phippery.collect

    :wd_location_column: A column specifying the unique
        integer location of the protein that
        the peptide is centered around

    :aa_sub_column: A column specifying the single aa char
        that we can expect to find in this peptide at the
        respective centered location

    :is_wt_column: A boolean column True/False specifying
        that the aa_sub at that loc on that protein is the
        wild type amino acid.
    """
    # TODO only for empirical samples

    ret = copy.deepcopy(ds)

    # differential selection is computed on a per "protein"/"loc" basis
    for protein, group_p in ds.groupby(ds.peptide_table.loc[:, protein_name_column]):
        for loc, group_p_l in group_p.groupby(
            group_p.peptide_table.loc[:, wd_location_column]
        ):

            wt_pep_id = group_p_l.peptide_id.where(
                group_p_l.peptide_table.loc[:, is_wt_column], drop=True
            ).peptide_id.values

            # sanity check 1, there should only be a single
            # wildtype codon per "Loc"/"Protein" combo
            assert len(wt_pep_id) == 1

            # diff selection is computed per sample.
            for sam_id in group_p_l.sample_id.values:

                # compute diff selection with each sample for all aa substitutions
                wt_enrichment = group_p_l.counts.loc[wt_pep_id[0], sam_id].values
                diff_sel = [
                    np.log2(e / wt_enrichment)
                    for e in group_p_l.counts.loc[:, sam_id].values
                ]
                ret.counts.loc[list(group_p_l.peptide_id.values), sam_id] = diff_sel

                # sanity check 2, the diff selection of the wildtype by def, is zero
                assert ret.counts.loc[wt_pep_id[0], sam_id] == 0.0

    return ret
