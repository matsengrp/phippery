"""
=================
Normalize
=================

A set of useful functions for normalizing enrichments
in a phippery dataset.
"""

import numpy as np
from numpy.linalg import svd
import xarray as xr
import pandas as pd
import itertools
import copy

from phippery.utils import iter_peptide_groups
from phippery.utils import iter_sample_groups
from phippery.utils import sample_id_coordinate_from_query
from phippery.utils import peptide_id_coordinate_from_query
from phippery.utils import get_annotation_table

# TODO J: example
def standardized_enrichment(
    ds,
    lib_ds,
    beads_ds,
    data_table="counts",
    inplace=True,
    new_table_name="std_enrichment"
):
    """Compute standardized enrichment of sample counts.
    This is the *fold enrichment* of each sample's frequency
    compared to the library average frequency, minus the mock IP
    (beads only control) frequency.

    Note
    ----
    Psuedo counts and exact calculation are
    derived from the bloom lab's formulated normalization
    heuristic for differential selection. See
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    lib_ds : xarray.DataSet
        The dataset of phage library control samples used in normalization

    beads_ds : xarray.DataSet
        The dataset of beads only control samples used in normalization

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray


    Returns
    -------
    
    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    for control in [lib_ds, beads_ds]:
        if type(control) != xr.Dataset:
            raise ValueError(
                "ds_lib_control_indicies must be of type list, even if there is only a single value"
            )

    std_enrichments = _comp_std_enr(
        counts=ds[data_table].to_pandas(),
        lib_counts=lib_ds[data_table].to_pandas(),
        mock_ip_counts=beads_ds[data_table].to_pandas(),
    )

    if inplace:
        ds[new_table_name] = xr.DataArray(std_enrichments)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(std_enrichments)
        return ds_copy


def _comp_std_enr(counts, lib_counts, mock_ip_counts):
    """Computes standardized enrichment."""

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


def enrichment(
    ds, lib_ds, data_table="counts", inplace=True, new_table_name="enrichment",
):
    """This function computes fold enrichment in the same fashion as
    the **standardized_enrichment**, but does *not* subtract beads only controls

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    lib_ds : xarray.DataSet
        The dataset of phage library control samples used in normalization

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
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
    """Compute enrichment values."""

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

# TODO J: add math
def svd_rank_reduction(
    ds,
    rank=1,
    data_table="enrichment",
    inplace=True,
    new_table_name="svd_rr",
):
    """
    This function computes the singular value decomposition, 
    then recomputes the enrichment matrix up to the rank specified.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    rank : int
        The number of ranks to include in the reconstruction.

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    low_rank_dt = copy.deepcopy(ds[data_table].to_pandas())

    # compute rank reduction decomposition matrices
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


# TODO J: col or column or feature, maybe?
def svd_aa_loc(
    ds,
    rank=1,
    data_table="enrichment",
    protein_name_column="Protein",
    location_col="Loc",
    aa_sub_col="aa_sub",
    inplace=True,
    new_table_name="svd_rr",
):
    """
    Compute singular value decomposition rank reduction
    on the aa / loc matrix by pivoting before computing decomposition
    and re-shaping to add to the dataset.

    Note
    ----

    This function is meant to be used specifically with phage-dms data
    where the peptide table includes a "loc" column, and an "aa_sub_col"
    which specifies the amino acid at that location.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    protein_name_column : str
        The peptide table feature which specifies which protein a specific peptide
        derives from.

    location_col : str
        The peptide table feature which specifies the site that a particular peptide
        is centered at.

    aa_sub_col : str
        The peptide table feature which specifies the amino acid at a given site.


    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    
    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
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

# TODO J: example
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
    """A generalized function to compute differential selection
    of amino acid variants in relation to the wildtype sequence.
    The function computed log fold change between enrichments
    of a wildtype and mutation at any given site.

    Note
    ----

    This function is meant to be used specifically with phage-dms data
    where the peptide table includes a "loc" column, and an "aa_sub_col"
    which specifies the amino acid at that location.

    Note
    ----
    This calculation of differential selection 
    is derived from the bloom lab's formulated from
    https://jbloomlab.github.io/dms_tools2/diffsel.html#id5


    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    scaled_by_wt : bool
        A boolean flag indicating whether or not you would like to multiply the
        differential selection value by

    smoothing_flank_size : int
        This parameter should only be used if **scaled_by_wt** is true.
        By specifying an integer greater than 0, you are scaling the differential
        selection value by enrichment values of the wildtype peptides surrounding
        , in both directions, a given site. The integer specified here then
        determines how many peptides are used for the scaling in both directions.

    groupby: list[str]
        This will specify which peptide feature groups such that site-mutation
        combinations are unique.

    loc_column : str
        The peptide table feature which specifies the site that a particular peptide
        is centered at. 

    is_wt_column : str
        The column specifying which peptides are wildtype.

    relu_bias : int
        If an integer is specified, then enrichment values less than
        1 are replaced by the specified value before computing differential
        selection.

    skip_samples : set
        sample id's which you do not want to calculate the differential selection on,
        such as controls. This function has many nested loops, so avoid computing
        on unnecessary samples.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray 

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    

    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    diff_sel = copy.deepcopy(ds[data_table])
    sw = smoothing_flank_size

    for group, group_ds in iter_peptide_groups(ds, groupby):

        wt_pep_id = peptide_id_coordinate_from_query(
            group_ds, [f"{is_wt_column} == True"]    
        )

        group_loc = group_ds.peptide_table.loc[wt_pep_id, loc_column].values
        for i, loc in enumerate(group_loc):

            loc_pid = peptide_id_coordinate_from_query(
                group_ds, [f"{loc_column} == {loc}"]
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
    Get a scalar from a wt sequence with a certain flank size.
    """

    if flank_size == 0:
        return wt_enr[i]

    lcase = i - flank_size < 0
    rcase = i + flank_size + 1 > len(wt_enr)
    lflank = wt_enr[i - flank_size : i] if not lcase else wt_enr[:i]
    rflank = wt_enr[i : i + flank_size + 1] if not rcase else wt_enr[i:]
    window_enr = lflank + rflank
    return sum(window_enr) / len(window_enr)


# TODO finish Doctring
def differential_selection_sample_groups(
    ds,
    sample_feature="library_batch",
    is_equal_to="batch_a",
    data_table="counts",
    aggregate_function=np.mean,
    inplace=True,
    new_table_name="sample_group_differential_selection",
):
    """This function computes differential selection
    between groups of samples.

    Note
    ----
    This function is still experimental.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    

    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

    diff_sel = copy.deepcopy(ds[data_table])
    group_id = sample_id_coordinate_from_query(ds, [f"{sample_feature} == '{is_equal_to}'"])
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
    a private function to compute differential selection of one values to a list of other values. 
    Optionally, you can scale each of the differential selection values by the base if desired.
    """

    if np.any(np.array(all_other_values) == 0):
        raise ZeroDivisionError(
            f"All values for which we are computing differential selection must be non-zero"
        )
    diff_sel = np.array([np.log2(v / base) for v in all_other_values])
    return diff_sel * scalar


def size_factors(
    ds, 
    inplace=True, 
    data_table="counts", 
    new_table_name="size_factors"
):
    """Compute size factors from Anders and Huber 2010.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    """

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


# TODO J: Add Math.
def counts_per_million(
        ds, 
        inplace=True, 
        new_table_name="cpm", 
        per_sample=True, 
        data_table="counts"
):
    """Compute counts per million.
    

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    per_sample : bool
        If True, compute counts per million separately for each sample.
        Otherwise, frequencies are computed as a ratio of each count to the sum of
        all counts in the dataset.

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 
    

    """

    if data_table not in ds:
        avail = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        raise KeyError(
            f"{data_table} is not included in dataset. \n available datasets: {avail}"
        )

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


# TODO J: Example, ascending vs decending? Should we do that?
def rank_data(
    ds, data_table="counts", inplace=True, per_sample=False, new_table_name=f"rank",
):
    """Compute the rank of specified enrichment layer.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to

    per_sample : bool
        If True, compute rank separately for each sample.
        Otherwise, frequencies are computed as a ratio of each count to the sum of
        all counts in the dataset.

    data_table : str
        The name of the enrichment layer you would like to fit mlxp to.

    new_table_name : str
        The name of the new layer you would like to append to the dataset.

    inplace : bool
        If True, then this function
        appends a dataArray to ds which is indexed with the same coordinate dimensions as
        'data_table'. If False, a copy of ds is returned with the appended dataArray
    

    Returns
    -------

    xarray.DataSet :
        If inplace is False, return a new DataSet object which has
        the enrichment values appended 

    """

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
