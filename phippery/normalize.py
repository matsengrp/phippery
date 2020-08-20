"""
@File: normalize.py

@Author: Jared Galloway

this file will contain functions which take
an xarray phip dataset and return a copy
where the counts have been normalized or transformed.
"""

import numpy as np
import xarray as xr
import pandas as pd
import itertools
import copy


def standardized_enrichment(
    ds,
    ds_lib_control_indices,
    ds_bead_control_indices,
    inplace=True,
    new_table_name="std_enrichment",
):
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
    # ret = copy.deepcopy(ds)
    std_enrichments = copy.deepcopy(ds.counts.to_pandas())

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
        std_enrichments.loc[:, sample_id] = sample_enrichment - pseudo_bead_enrichment

    if inplace:
        ds[new_table_name] = xr.DataArray(std_enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(std_enrichments)
        return ds_copy


def enrichment(ds, ds_lib_control_indices, inplace=True, new_table_name="enrichment"):
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

    # we are going to add an augmented counts matrix
    enrichments = copy.deepcopy(ds.counts.to_pandas())

    # find controls and average all
    ds_lib_counts_mean = enrichments[ds_lib_control_indices].mean(axis=1)
    ds_lib_counts_mean_sum = sum(ds_lib_counts_mean)

    # compute all sample standardized enrichment
    for sample_id, sample in enrichments.iteritems():
        if sample_id in list(ds_lib_control_indices):
            continue
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_control = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        enrichments.loc[:, sample_id] = sample_enrichment

    if inplace:
        ds[new_table_name] = xr.DataArray(enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(enrichments)
        return ds_copy


def differential_selection(
    ds,
    scaled_by_wt=False,
    protein_name_column="Protein",
    wd_location_column="Loc",
    aa_sub_column="aa_sub",
    is_wt_column="is_wt",
    inplace=True,
    new_table_name="differential_selection",
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

    # TODO actually, we could in principal do this
    # on any of the data tables.
    if "enrichment" not in ds:
        raise KeyError(
            "enrichment matrix must be computed before differential expression"
        )
    # ret = copy.deepcopy(ds)
    # TODO use numpy array_like
    diff_sel = copy.deepcopy(ds.enrichment)

    # TODO use iter method from utils
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

                # TODO this should either be a parameter or
                # we should warn heavily against putting controls in?
                # maybe this doesn't matter
                emp = ds.sample_table.loc[sam_id, "control_status"].values
                if emp != "empirical":
                    continue

                # compute diff selection with each sample for all aa substitutions
                wt_enrichment = group_p_l.enrichment.loc[wt_pep_id[0], sam_id].values
                local_diff_sel = [
                    np.log2(e / wt_enrichment)
                    for e in group_p_l.enrichment.loc[:, sam_id].values
                ]
                scaled = (
                    local_diff_sel
                    if not scaled_by_wt
                    else local_diff_sel * wt_enrichment
                )
                # ret.counts.loc[list(group_p_l.peptide_id.values), sam_id] = scaled
                diff_sel.loc[list(group_p_l.peptide_id.values), sam_id] = scaled

                # sanity check 2, the diff selection of the wildtype by def, is zero
                # assert ret.counts.loc[wt_pep_id[0], sam_id] == 0.0
                assert diff_sel.loc[wt_pep_id[0], sam_id] == 0.0

    if inplace:
        ds[new_table_name] = diff_sel
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = diff_sel
        return ds_copy


def size_factors(ds, inplace=True, new_table_name="size_factors"):
    """Compute size factors from Anders and Huber 2010

    counts is a numpy array

    This function was originally implemented in phip-stat here:
    https://github.com/lasersonlab/phip-stat
    """

    # TODO use numpy array_like
    size_factors = copy.deepcopy(ds.counts.to_pandas())
    counts = size_factors.values

    masked = np.ma.masked_equal(counts, 0)
    geom_means = (
        np.ma.exp(np.ma.log(masked).sum(axis=1) / (~masked.mask).sum(axis=1))
        .data[np.newaxis]
        .T
    )

    size_factors = (
        size_factors / np.ma.median(masked / geom_means, axis=0).data
    ).round(2)

    if inplace:
        ds[new_table_name] = xr.DataArray(size_factors)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(size_factors)
        return ds_copy


def cpm(ds, inplace=True, new_table_name="cpm"):
    """compute counts per million for the given data
    and then add it to the dataset as a new table"""

    # TODO use numpy array_like
    new = copy.deepcopy(ds.counts.to_pandas())
    cpm = (new / (new.sum() / 1e6)).round(2)

    if inplace:
        ds[new_table_name] = xr.DataArray(cpm)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(cpm)
        return ds_copy


def rank_data_per_sample(
    ds, data_table="counts", inplace=True, new_table_name=f"sample_rank"
):
    """given a data set and a table,
    compute the rank of each sample's peptide
    score wrt the data_table. Add this rank table
    to the dataset"""

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    # TODO: use array_like we know hthe dt should be int
    new = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    for sid in ds.sample_id.values:
        sample_data = ds[f"{data_table}"].loc[:, sid].values
        temp = sample_data.argsort()
        ranks = np.empty_like(temp)
        ranks[temp] = np.arange(len(sample_data))
        new.loc[:, sid] = ranks.flatten()

    if inplace:
        ds[new_table_name] = xr.DataArray(new)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(new)
        return ds_copy


def rank_data_table(ds, data_table, inplace=True, new_table_name=f"table_rank"):
    """given a data set and a table,
    compute the rank of every sample_peptide combination.
    score wrt the data_table. Add this rank table
    to the dataset"""

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    new = copy.deepcopy(ds[f"{data_table}"].to_pandas())
    sample_data_sh = ds[f"{data_table}"].values.shape
    sample_data = ds[f"{data_table}"].values.flatten()
    temp = sample_data.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(sample_data))
    new.loc[:, :] = ranks.reshape(sample_data_sh)

    if inplace:
        ds[new_table_name] = xr.DataArray(new)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(new)
        return ds_copy
