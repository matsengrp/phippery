"""
@File: normalize.py

@Author: Jared Galloway

this file will contain functions which take
an xarray phip dataset and return a copy
where the counts have been normalized or transformed.
"""

import numpy as np
from numpy.linalg import svd
import xarray as xr
import pandas as pd
import itertools
import copy

from phippery.utils import iter_peptide_groups
from phippery.utils import iter_sample_groups
from phippery.utils import id_coordinate_subset
from phippery.tidy import tidy_ds


def standardized_enrichment(
    ds,
    lib_ds,
    mock_ip_ds,
    data_table="counts",
    inplace=True,
    new_table_name="std_enrichment",
):
    """
    return a new xarray dataset same as the input
    except with the counts converted to standard enrichment.

    pseudo counts are added like so:
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    if inplace other values will not be tampered with.

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
    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    # TODO and verify(ds) == True
    for control in [lib_ds, mock_ip_ds]:
        if type(control) != xr.Dataset:
            raise ValueError(
                "ds_lib_control_indicies must be of type list, even if there is only a single value"
            )

    std_enrichments = _comp_std_enr(
        counts=ds[data_table].to_pandas(),
        lib_counts=lib_ds[data_table].to_pandas(),
        mock_ip_counts=mock_ip_ds[data_table].to_pandas(),
    )

    if inplace:
        ds[new_table_name] = xr.DataArray(std_enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(std_enrichments)
        return ds_copy


def _comp_std_enr(counts, lib_counts, mock_ip_counts):

    normalized_ds_counts = copy.deepcopy(counts)

    # find controls and average all
    ds_lib_counts_mean = lib_counts.mean(axis=1)
    ds_bead_counts_mean = mock_ip_counts.mean(axis=1)
    ds_lib_counts_mean_sum = sum(ds_lib_counts_mean)

    # compute beads control std enrichment
    pseudo_sample = ds_bead_counts_mean + max(
        1, sum(ds_bead_counts_mean) / ds_lib_counts_mean_sum
    )
    pseudo_lib_counts = ds_lib_counts_mean + max(
        1, ds_lib_counts_mean_sum / sum(ds_bead_counts_mean)
    )
    pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
    pseudo_lib_counts_freq = pseudo_lib_counts / sum(pseudo_lib_counts)
    pseudo_bead_enrichment = pseudo_sample_freq / pseudo_lib_counts_freq

    # compute all sample standardized enrichment
    for sample_id, sample in counts.iteritems():
        pseudo_sample = sample + max(1, sum(sample) / ds_lib_counts_mean_sum)
        pseudo_lib_counts = ds_lib_counts_mean + max(
            1, ds_lib_counts_mean_sum / sum(sample)
        )
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_counts_freq = pseudo_lib_counts / sum(pseudo_lib_counts)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_counts_freq
        normalized_ds_counts.loc[:, sample_id] = (
            sample_enrichment - pseudo_bead_enrichment
        )

    return normalized_ds_counts


# def standardized_enrichment(
#    ds,
#    lib_ds,
#    mock_ip_ds,
#    data_table="counts",
#    inplace=True,
#    new_table_name="std_enrichment",
# ):
def enrichment(
    ds, lib_ds, data_table="counts", inplace=True, new_table_name="enrichment",
):
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

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    if type(lib_ds) != xr.Dataset:
        raise ValueError(
            "ds_lib_control_indicies must be of type list, even if there is only a single value"
        )

    enrichments = _comp_enr(
        counts=ds[data_table].to_pandas(), lib_counts=lib_ds[data_table].to_pandas()
    )

    if inplace:
        ds[new_table_name] = xr.DataArray(enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(enrichments)
        return ds_copy


def _comp_enr(counts, lib_counts):
    # we are going to add an augmented counts matrix
    enrichments = copy.deepcopy(counts)

    # find controls and average all
    lib_counts_mean = lib_counts.mean(axis=1)
    lib_counts_mean_sum = sum(lib_counts_mean)

    # compute all sample standardized enrichment
    for sample_id, sample in enrichments.iteritems():

        pseudo_sample = sample + max(1, sum(sample) / lib_counts_mean_sum)
        pseudo_lib_control = lib_counts_mean + max(1, lib_counts_mean_sum / sum(sample))
        pseudo_sample_freq = pseudo_sample / sum(pseudo_sample)
        pseudo_lib_control_freq = pseudo_lib_control / sum(pseudo_lib_control)
        sample_enrichment = pseudo_sample_freq / pseudo_lib_control_freq
        enrichments.loc[:, sample_id] = sample_enrichment

    return enrichments


def svd_rank_reduction(
    ds,
    rank=1,
    data_table="enrichment",
    scaled_by_wt=False,
    inplace=True,
    new_table_name="svd_rr",
):
    """
    compute singular value decomposition rank reduction
    on any data table in the dataset. Add the resulting
    rank reduced layer to the dataset.

    :param: r <int> Number of ranks in re-composition estimate.
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    low_rank_dt = copy.deepcopy(ds[data_table].to_pandas())

    # compute rank reduction decompisition matrices
    U, S, V = svd(low_rank_dt.values)

    # Grab the first X outer products in the finite summation of rank layers.
    low_rank = U[:, :rank] @ np.diag(S[:rank]) @ V[:rank, :]
    low_rank_dt.loc[:, :] = low_rank
    svd_rr_approx = xr.DataArray(low_rank_dt, dims=ds[data_table].dims)

    if inplace:
        ds[new_table_name] = svd_rr_approx
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = svd_rr_approx
        return ds_copy


def svd_aa_loc(
    ds,
    rank=1,
    data_table="enrichment",
    scaled_by_wt=False,
    protein_name_column="Protein",
    location_col="Loc",
    aa_sub_col="aa_sub",
    inplace=True,
    new_table_name="svd_rr",
):
    """
    compute singular value decomposition rank reduction
    on the aa / loc matrix by pivoting before computing decomposiion
    and re-shaping to add to the dataset.

    :param: r <int> Number of ranks in re-composition estimate.
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    low_rank_dt = copy.deepcopy(ds[data_table].to_pandas())

    for prot, prot_ds in iter_peptide_groups(ds, protein_name_column):

        for sid in prot_ds.sample_id.values:

            # grab the single sample ds
            rep_ds = prot_ds.loc[
                dict(
                    sample_id=[sid],
                    sample_metadata=["sample_ID"],
                    peptide_metadata=[aa_sub_col, location_col],
                )
            ]

            # melt
            tidy = tidy_ds(rep_ds)

            # Pivot so that we get the (aa X Loc)
            piv = tidy.pivot(index=aa_sub_col, columns=location_col, values=data_table)

            # Preserve the indices for population of new enrichment table
            piv_index = tidy.pivot(
                index=aa_sub_col, columns=location_col, values="peptide_id"
            )

            # compute rank reduction decompisition matrices
            U, S, V = svd(piv)

            # Grab the first X outer products in the finite summation of rank layers.
            low_rank = U[:, :rank] @ np.diag(S[:rank]) @ V[:rank, :]

            low_rank_piv = pd.DataFrame(low_rank, index=piv.index, columns=piv.columns)
            melted_values = pd.melt(low_rank_piv.reset_index(), id_vars=[aa_sub_col])
            melted_index = pd.melt(piv_index.reset_index(), id_vars=[aa_sub_col])
            melted_values["peptide_id"] = melted_index["value"]
            low_rank_dt.loc[melted_values["peptide_id"], sid] = melted_values[
                "value"
            ].values

    svd_rr_approx = xr.DataArray(low_rank_dt, dims=ds[data_table].dims)

    if inplace:
        ds[new_table_name] = svd_rr_approx
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = svd_rr_approx
        return ds_copy


def differential_selection_wt_mut(
    ds,
    data_table="enrichment",
    scaled_by_wt=False,
    smoothing_flank_size=0,
    groupby=["Protein"],
    loc_column="Loc",
    is_wt_column="is_wt",
    inplace=True,
    new_table_name="wt_mutant_differential_selection",
    relu_bias=None,
    skip_samples=set(),
):
    # TODO
    """
    A generalized function to compute differential selection
    of amino acid variants in relation to the wildtype sequence.
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    diff_sel = copy.deepcopy(ds[data_table])
    sw = smoothing_flank_size

    # iterate through groups which have a unique loc
    # from tqdm import tqdm
    #for group, group_ds in tqdm(iter_peptide_groups(ds, groupby)):
    for group, group_ds in iter_peptide_groups(ds, groupby):

        wt_pep_id = id_coordinate_subset(
            group_ds, table="peptide_table", where=is_wt_column, is_equal_to=True,
        )

        group_loc = group_ds.peptide_table.loc[wt_pep_id, loc_column].values
        #for i, loc in tqdm(enumerate(group_loc), leave=False):
        for i, loc in enumerate(group_loc):

            loc_pid = id_coordinate_subset(
                group_ds, table="peptide_table", where=loc_column, is_equal_to=loc
            )
            loc_ds = group_ds.loc[dict(peptide_id=loc_pid)]

            # check that skip samples is of type list
            sams = set(loc_ds.sample_id.values) - set(skip_samples)
            for sam_id in sams:

                wt_seq_enr = group_ds[data_table].loc[wt_pep_id, sam_id].values
                wt_enrichment = float(wt_seq_enr[i])
                scalar = (
                    _wt_window_scalar(list(wt_seq_enr), i, sw) if scaled_by_wt else 1
                )
                values = loc_ds[data_table].loc[:, sam_id].values

                if relu_bias is not None:
                    values[values < 1] = relu_bias
                dsel = _comp_diff_sel(wt_enrichment, values, scalar)

                diff_sel.loc[list(loc_ds.peptide_id.values), sam_id] = dsel
                assert diff_sel.loc[wt_pep_id[0], sam_id] == 0.0

    if inplace:
        ds[new_table_name] = diff_sel
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = diff_sel
        return ds_copy


def _wt_window_scalar(wt_enr, i, flank_size):
    """
    get a scalar from a wt sequence with a certain flank size.
    """

    if flank_size == 0:
        return wt_enr[i]

    lcase = i - flank_size < 0
    rcase = i + flank_size + 1 > len(wt_enr)
    lflank = wt_enr[i - flank_size : i] if not lcase else wt_enr[:i]
    rflank = wt_enr[i : i + flank_size + 1] if not rcase else wt_enr[i:]
    window_enr = lflank + rflank
    return sum(window_enr) / len(window_enr)


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
    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

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


def _comp_diff_sel(base, all_other_values, scalar=1):
    """
    a private function to compute differential selection of one values to a list of other values. Optionally, you can scale each of the differential selection values by the base if desired.
    """

    if np.any(np.array(all_other_values) == 0):
        raise ZeroDivisionError(
            f"All values for which we are computing differential selection must be non-zero"
        )
    diff_sel = np.array([np.log2(v / base) for v in all_other_values])
    # return diff_sel if not scaled_by_base else diff_sel * base
    return diff_sel * scalar


def size_factors(ds, inplace=True, data_table="counts", new_table_name="size_factors"):
    """Compute size factors from Anders and Huber 2010"""

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

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

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

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
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

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
