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

from phippery.utils import iter_peptide_groups
from phippery.utils import iter_sample_groups
from phippery.utils import id_coordinate_subset


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

    if type(ds_lib_control_indices) != list:
        raise ValueError(
            "ds_lib_control_indicies must be of type list, even if there is only a single values"
        )

    if type(ds_bead_control_indices) != list:
        raise ValueError(
            "ds_bead_control_indicies must be of type list, even if there is only a single values"
        )

    ds_counts = ds.counts.to_pandas()
    std_enrichments = _comp_std_enr(
        ds_counts, ds_lib_control_indices, ds_bead_control_indices
    )

    if inplace:
        ds[new_table_name] = xr.DataArray(std_enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(std_enrichments)
        return ds_copy


def _comp_std_enr(counts, lib_controls, beads_controls):

    std_enrichments = copy.deepcopy(counts)

    # find controls and average all
    ds_lib_counts_mean = counts[lib_controls].mean(axis=1)
    ds_bead_counts_mean = counts[beads_controls].mean(axis=1)
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
    for sample_id, sample in counts.iteritems():
        if sample_id in list(lib_controls) + list(beads_controls):
            continue
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_control = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        std_enrichments.loc[:, sample_id] = sample_enrichment - pseudo_bead_enrichment

    return std_enrichments


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
    if type(ds_lib_control_indices) != list:
        raise ValueError(
            "ds_lib_control_indicies must be of type list, even if there is only a single values"
        )

    ds_counts = ds.counts.to_pandas()
    enrichments = _comp_enr(ds_counts, ds_lib_control_indices,)

    if inplace:
        ds[new_table_name] = xr.DataArray(enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(enrichments)
        return ds_copy


def _comp_enr(counts_df, lib_controls):
    # we are going to add an augmented counts matrix
    enrichments = copy.deepcopy(counts_df)

    # find controls and average all
    ds_lib_counts_mean = enrichments[lib_controls].mean(axis=1)
    ds_lib_counts_mean_sum = sum(ds_lib_counts_mean)

    # compute all sample standardized enrichment
    for sample_id, sample in enrichments.iteritems():
        if sample_id in list(lib_controls):
            continue
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_control = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        enrichments.loc[:, sample_id] = sample_enrichment

    return enrichments


def differential_selection_wt_mut(
    ds,
    data_table="enrichment",
    scaled_by_wt=False,
    protein_name_column="Protein",
    wd_location_column="Loc",
    is_wt_column="is_wt",
    inplace=True,
    new_table_name="wt_mutant_differential_selection",
):

    diff_sel = copy.deepcopy(ds[data_table])
    for protein, protein_ds in iter_peptide_groups(ds, protein_name_column):
        for loc, protein_loc_ds in iter_peptide_groups(protein_ds, wd_location_column):
            wt_pep_id = id_coordinate_subset(
                protein_loc_ds,
                table="peptide_table",
                where=is_wt_column,
                is_equal_to=True,
            )
            assert len(wt_pep_id) == 1

            for sam_id in protein_loc_ds.sample_id.values:

                wt_enrichment = (
                    protein_loc_ds[data_table].loc[wt_pep_id[0], sam_id].values
                )
                values = protein_loc_ds[data_table].loc[:, sam_id].values
                dsel = _comp_diff_sel(wt_enrichment, values, scaled_by_wt)

                diff_sel.loc[list(protein_loc_ds.peptide_id.values), sam_id] = dsel
                assert diff_sel.loc[wt_pep_id[0], sam_id] == 0.0

    if inplace:
        ds[new_table_name] = diff_sel
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = diff_sel
        return ds_copy


def differential_selection_sample_groups(
    ds,
    sample_feature="library_batch",
    is_equal_to="batch_a",
    data_table="counts",
    aggregate_function=np.mean,
    inplace=True,
    new_table_name="sample_group_differential_selection",
):
    """
    This function should compute differential selection
    between groups of samples rather than wildtype vs mutant
    """

    # TODO Write Checks here

    diff_sel = copy.deepcopy(ds[data_table])
    group_id = id_coordinate_subset(ds, where=sample_feature, is_equal_to=is_equal_to)
    group_enrichments = ds[data_table].loc[:, group_id].values
    group_agg = np.apply_along_axis(aggregate_function, 1, group_enrichments)
    for agg_enrich, peptide_id in zip(group_agg, ds.peptide_id.values):
        all_other_sample_values = ds[data_table].loc[peptide_id, :].values
        peptide_diff_sel = _comp_diff_sel(agg_enrich, all_other_sample_values)
        diff_sel.loc[peptide_id, :] = peptide_diff_sel

    if inplace:
        ds[new_table_name] = diff_sel
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = diff_sel
        return ds_copy


def _comp_diff_sel(base, all_other_values, scaled_by_base=False):
    """
    a private function to compute differential selection of one values to a list of other values. Optionally, you can scale each of the differential selection values by the base if desired.
    """

    if np.any(np.array(all_other_values) == 0):
        raise ZeroDivisionError(
            f"All values for which we are computing differential selection must be non-zero"
        )
    diff_sel = np.array([np.log2(v / base) for v in all_other_values])
    return diff_sel if not scaled_by_base else diff_sel * base


def size_factors(ds, inplace=True, data_table="counts", new_table_name="size_factors"):
    """Compute size factors from Anders and Huber 2010"""

    size_factors = _comp_size_factors(ds[data_table].to_pandas().values)
    sf_da = xr.DataArray(size_factors, dims=ds[data_table].dims)

    if inplace:
        ds[new_table_name] = sf_da
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = sf_da
        return ds_copy


def _comp_size_factors(counts):

    size_factors = copy.deepcopy(counts)
    masked = np.ma.masked_equal(counts, 0)

    if type(masked.mask) != np.ndarray:
        bool_mask = np.full(counts.shape, False, dtype=bool)
    else:
        bool_mask = ~masked.mask

    geom_means = (
        np.ma.exp(np.ma.log(masked).sum(axis=1) / (bool_mask).sum(axis=1))
        .data[np.newaxis]
        .T
    )

    size_factors = (
        size_factors / np.ma.median(masked / geom_means, axis=0).data
    ).round(2)

    return size_factors


def cpm(ds, inplace=True, new_table_name="cpm", per_sample=False, data_table="counts"):
    """compute counts per million for the given data
    and then add it to the dataset as a new table"""

    # TODO use numpy array_like
    counts = ds[data_table].to_pandas().values
    cpm = _comp_cpm_per_sample(counts) if per_sample else _comp_cpm(counts)
    cpm_da = xr.DataArray(cpm, dims=ds[data_table].dims)

    if inplace:
        ds[new_table_name] = cpm_da
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = cpm_da
        return ds_copy


def _comp_cpm(counts):

    ret = copy.deepcopy(counts)
    return (ret / (ret.sum() / 1e6)).round(2)


def _comp_cpm_per_sample(counts):

    ret = copy.deepcopy(counts)
    return (ret / (ret.sum(axis=0) / 1e6)).round(2)


def rank_data(
    ds, data_table="counts", inplace=True, per_sample=False, new_table_name=f"rank",
):
    """given a data set and a table,
    compute the rank of each sample's peptide
    score wrt the data_table. Add this rank table
    to the dataset"""

    if data_table not in ds:
        raise KeyError(f"{data_table} is not included in dataset.")

    counts = ds[data_table].to_pandas().values
    cpm = _comp_rank_per_sample(counts) if per_sample else _comp_rank(counts)
    cpm_da = xr.DataArray(cpm, dims=ds[data_table].dims)

    if inplace:
        ds[new_table_name] = cpm_da
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = cpm_da
        return ds_copy


def _comp_rank(enrichment):

    ret = np.ones(enrichment.shape)
    sample_data_sh = enrichment.shape
    sample_data = enrichment.flatten()
    temp = sample_data.argsort()
    ranks = np.empty_like(temp)
    ranks[temp] = np.arange(len(sample_data))
    ret[:, :] = ranks.reshape(sample_data_sh)

    return ret.astype(int)


def _comp_rank_per_sample(enrichment):

    ret = np.ones(enrichment.shape)
    for sid in range(enrichment.shape[1]):
        sample_data = enrichment[:, sid]
        temp = sample_data.argsort()
        ranks = np.empty_like(temp)
        ranks[temp] = np.arange(len(sample_data))
        ret[:, sid] = ranks.flatten()

    return ret.astype(int)
