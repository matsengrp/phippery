"""
@File: utils.py

@Author: Jared Galloway

This file will include some helpful functions
for the phippery package CLI. The primary
data struct we put counts data into is
"""

# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as st

# built-in python3
import os
import re
import copy
import itertools
import pickle
from functools import reduce
from collections import defaultdict


def get_all_sample_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a sample table column
    """

    all_exp = ds.sample_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def get_all_peptide_metadata_factors(ds, feature):
    """
    return a list of all available factors in
    a peptide table column
    """

    all_exp = ds.peptide_table.loc[:, feature]
    return [x for x in set(all_exp.values) if x == x]


def convert_peptide_table_to_fasta(peptide_table, out):
    """
    Take in peptide metadata dataframe, and write a fasta
    format representation of the oligos
    """

    fasta_fp = open(out, "w")
    peptide_table = pd.read_csv(peptide_table, index_col=0, header=0)
    requirements = ["Oligo"]
    assert peptide_table.index.name == "peptide_id"
    assert np.all([x in peptide_table.columns for x in requirements])
    for index, row in peptide_table.iterrows():
        ref_sequence = trim_index(row["Oligo"])
        fasta_fp.write(f">{index}\n{ref_sequence}\n")


def trim_index(sequence):
    return "".join([nt for nt in sequence if nt.isupper()])


#######################################################
# SLICING/INDEX FUNCTIONS
#######################################################


def iter_groups(ds, by, dim="sample"):
    """
    returns an iterator yeilding subsets of the provided dataset,
    grouped by items in the metadata of either dimension.
    note that this provides references to subsets of the data,
    and thus modifying the yeilded objects directly modifies the
    original dataset being passed in.
    """

    table = get_annotation_table(ds, dim=dim)
    for group, group_df in table.groupby(by):
        group_ds = ds.loc[{f"{dim}_id": list(group_df.index.values)}]
        yield group, group_ds


def iter_sample_groups(*args):
    """
    DEPRECATED - use 'iter_groups()' instead

    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the sample metadata coodinate.
    """

    return iter_groups(*args, "sample")


def iter_peptide_groups(*args):

    """DEPRECATED - use 'iter_groups()' instead

    returns an iterator yeilding subsets of the provided dataset,
    grouped by an item on the peptide metadata coodinate. """

    return iter_groups(*args, "peptide")


def id_coordinate_from_query(ds, query_df):

    """
    Given a df with columns 'dimension' and 
    """

    # st = ds.sample_table.to_pandas().infer_objects()
    sq = list(query_df.loc[query_df["dimension"] == "sample", "expression"].values)
    sid = sample_id_coordinate_from_query(ds, sq)

    # pt = ds.peptide_table.to_pandas().infer_objects()
    pq = list(query_df.loc[query_df["dimension"] == "peptide", "expression"].values)
    pid = peptide_id_coordinate_from_query(ds, pq)

    return sid, pid


def peptide_id_coordinate_from_query(ds, query_list: list, *args, **kwargs):
    """Take in a list of queries and return the peptide id index resulting
    from query """

    if len(query_list) == 0:
        return list(ds.peptide_id.values)

    peptide_table = get_peptide_table(ds)
    return list(peptide_table.query(" & ".join(query_list)).index.values)


def sample_id_coordinate_from_query(ds, query_list: list, *args, **kwargs):
    """Take in a list of queries and return the sample id index resulting
    from query """

    if len(query_list) == 0:
        return list(ds.sample_id.values)

    sample_table = get_sample_table(ds)
    return list(sample_table.query(" & ".join(query_list)).index.values)


#######################################################
# XARRAY MANIPULATION
#######################################################


def load(path):
    """
    simple wrapper for loading xarray datasets from pickle binary
    """

    return pickle.load(open(path, "rb"))


def dump(ds, path):
    """
    simple wrapper for dump'ing xarray datasets to pickle binary
    """

    pickle.dump(ds, open(path, "wb"))


def get_annotation_table(ds, dim="sample"):
    """return a copy of the peptide table after converting all
    the datatypes applying the pandas NaN hueristic"""

    st = (ds[f"{dim}_table"]
        .to_pandas()
        .convert_dtypes()
    )
    return st


def get_sample_table(*args):
    """DEPRECATED - use 'get_annotation_table()' instead

    return a copy of the peptide table after converting all
    the datatypes using the pandas NaN hueristic"""

    return get_annotation_table(*args, dim = "sample")


def get_peptide_table(*args):
    """DEPRECATED - use 'get_annotation_table()' instead
    
    return a copy of the sample table after converting all
    the datatypes using the pandas NaN heuristic"""

    return get_annotation_table(*args, dim = "peptide")


def stitch_dataset(
    counts, peptide_table, sample_table,
):
    """Simply stitch together the xarray matrix
    given the enrichments, sample metadata, and peptide metadata
    """

    # make sure the coordinated match up.
    assert np.all(counts.columns == sample_table.index)
    assert np.all(counts.index == peptide_table.index)

    # we are returning the xarray dataset organized by four coordinates seen below.
    pds = xr.Dataset(
        {
            "counts": (["peptide_id", "sample_id"], counts),
            "sample_table": (["sample_id", "sample_metadata"], sample_table),
            "peptide_table": (["peptide_id", "peptide_metadata"], peptide_table),
        },
        coords={
            "sample_id": counts.columns.values,
            "peptide_id": counts.index.values,
            "sample_metadata": sample_table.columns.values,
            "peptide_metadata": peptide_table.columns.values,
        },
    )
    return pds


def collect_merge_prune_count_data(counts):
    """ DEPRECATED - please use collect_counts() instead """

    return collect_counts(counts)


def collect_counts(counts):
    """
    This function takes in a list of paths which
    contains the counts for each peptide alignment
    for each sample. These files should contain
    no header.

    :param: counts <str> - a list of paths leading
    to raw peptide enrichment counts for each sample
    """

    load = lambda path, sample: pd.read_csv(  # noqa
        path, index_col=0, sep="\t", names=["peptide_id", sample]
    )

    sample_dataframes = [
        load(path, int(os.path.basename(path).split(".")[0])) for path in counts
    ]

    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        sample_dataframes,
    ).fillna(0)

    merged_counts_df.columns = merged_counts_df.columns.astype(int)
    merged_counts_df.index = merged_counts_df.index.astype(int)
    merged_counts_df.sort_index(inplace=True)
    merged_counts_df.sort_index(axis=1, inplace=True)

    return merged_counts_df


def tidy_ds(*args, **kwargs):
    """DEPRECATED - please use to_tall_csv()"""
    return to_tall(*args, **kwargs)


def to_tall(ds):
    """
    return the provided dataset in 'tall form'.

    This means that all information from the entire dataset and
    all the tables it includes will be 'melted' into a single dataframe.
    This format is far less effecient in terms of storage, and should not
    be used on the entire dataset, but rather, some subset of the dataset
    when you are ready to plot the relevent information.

    The means that for each row, you will get the following columns:

    sample_id,
    peptide_id,
    all sample metadata (a column for each),
    all peptide metadata (a column for each),
    all enrichment tables (a column for each).

    Ideal for gg plotting.
    """

    # melt all data tables in the dataset
    melted_data = [
        ds[f"{dt}"]
        .to_pandas()
        .reset_index()
        .melt(id_vars=["peptide_id"])
        .rename({"value": f"{dt}"}, axis=1)
        for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    ]

    # merge all melted tables so each get's a column in the final df
    # sample_coord_dim = ds.attrs["sample_coord_dim"]
    merged_counts_df = reduce(
        lambda l, r: pd.merge(l, r, on=["peptide_id", "sample_id"]), melted_data
    )

    # grab only the columns of the metadata we want in the resulting dataframe
    peptide_table = ds.peptide_table.to_pandas().reset_index().infer_objects()
    sample_table = ds.sample_table.to_pandas().reset_index().infer_objects()

    # merge the metadata into the melted datatables
    data_peptide = merged_counts_df.merge(peptide_table, on="peptide_id")
    return data_peptide.merge(sample_table, on="sample_id")


def to_wide(ds, file_prefix):
    """write an xarray object to wide style csv's
    with a specific prefix
    """

    ret = {}
    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    for dt in enr_layers:
        layer = copy.deepcopy(ds[f"{dt}"].to_pandas())
        layer.index.name = ""
        ret[dt] = layer

    for at in ["sample", "peptide"]:
        ret[at] = get_annotation_table(ds, dim=at)

    return ret


def dataset_to_wide_csv(*args, **kwargs):
    """DEPRECATED - please use to_wide_csv()"""
    return to_wide_csv(*args, **kwargs)


def to_wide_csv(ds, file_prefix):
    """write an xarray object to wide style csv's
    with a specific prefix
    """

    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    for dt in enr_layers:
        layer = ds[f"{dt}"].to_pandas()
        layer.index.name = ""
        layer.to_csv(
                f"{file_prefix}_{dt}.csv", 
                na_rep="NA", 
                index_label=None
        )

    for at in ["sample", "peptide"]:
        get_annotation_table(ds, dim=at).to_csv(
                f"{file_prefix}_{at}_annotation_table.csv"
        )


def dataset_from_csv(
        counts_matrix_filename, 
        peptide_table_filename, 
        sample_table_filename
):
    """
    """

    counts_df = collect_counts_matrix(counts_matrix_filename)
    peptide_df = collect_peptide_table(peptide_table_filename)
    sample_df = collect_sample_table(sample_table_filename)

    return stitch_dataset(
        counts = counts_df, 
        peptide_table=peptide_df, 
        sample_table=sample_df,
    )


def collect_sample_table(sample_table_filename: str):
    """
    """

    sample_table = pd.read_csv(
            sample_table_filename, 
            sep=",",
            index_col=0, 
            header=0
    ) 


    if sample_table.index.name != 'sample_id':
        raise ValueError("The name of the index must be 'sample_id'")

    if sample_table.index.dtype != 'int64':
        raise ValueError("The index values for sample_id must be inferred as integers")

    sample_table.sort_index(inplace=True)
    return sample_table


def collect_peptide_table(peptide_table_filename: str):
    """
    """

    peptide_table = pd.read_csv(
            peptide_table_filename, 
            sep=",", 
            index_col=0, 
            header=0
    )

    if peptide_table.index.name != 'peptide_id':
        raise ValueError

    if peptide_table.index.dtype != 'int64':
        raise ValueError("The index values for peptide_id must be inferred as integers")

    peptide_table.sort_index(inplace=True)
    return peptide_table


def collect_counts_matrix(counts_matrix_filename: str):
    """
    """

    counts_matrix = pd.read_csv(
        counts_matrix_filename, 
        sep=",", 
        index_col=0, 
        header=0
    )

    try:
        counts_matrix.columns = counts_matrix.columns.astype(int)
    except:
        raise ValueError("column header values much be able to cast to type 'int' to match peptide table index")

    try:
        counts_matrix.index = counts_matrix.index.astype(int)
    except:
        raise ValueError("row index values much be able to cast to type 'int' to match peptide table index")

    counts_matrix.sort_index(inplace=True)
    counts_matrix.sort_index(axis=1, inplace=True)

    return counts_matrix


#######################################################
# Collapse Functionality
#######################################################


def throw_mismatched_features(df, by):

    """When you collapse by some set of columns in the dataframe,
    keep only features which homogeneous within groups.

    This is similar to 'DataFrameGroupby.first()' method,
    but instead of keeping the first factor level appearing for each group
    category, we only throw any features which are note homogeneous within groups.
    You could achieve the same functionality by listing features you know to be
    homogeneous in the 'by' parameter."""

    # Find out which collapse features are shared within groups
    collapsed_sample_metadata = defaultdict(list)
    for i, (group, group_df) in enumerate(df.groupby(by)):
        for column, value in group_df.iteritems():
            v = value.values
            if np.all(v == v[0]) or np.all([n != n for n in v]):
                collapsed_sample_metadata[column].append(v[0])

    # Throw out features that are not shared between groups
    to_throw = [
        key for key, value in collapsed_sample_metadata.items() if len(value) < i + 1
    ]
    [collapsed_sample_metadata.pop(key) for key in to_throw]
    return pd.DataFrame(collapsed_sample_metadata)


def mean_pw_cc_by_multiple_tables(ds, by, dim="sample", data_tables="all"):

    """A wrapper for computing pw cc within groups defined
    with the 'by' parameter. Merges the correlations into a
    single table"""

    # Compute mean pw cc on all possible data tables
    if data_tables == "all":
        data_tables = list(set(ds.data_vars) - set(["sample_table", "peptide_table"]))

    # Some error handling
    if dim not in ["sample", "peptide"]:
        raise ValueError(f"parameter 'dim' must be either 'sample' or 'peptide'")

    groups_avail = ds[f"{dim}_metadata"].values
    for data_table in data_tables:
        if data_table not in ds:
            raise KeyError(f"{data_table} is not included in dataset.")

    for group in by:
        if group not in groups_avail:
            raise KeyError(
                f"{group} is not included as a column in the {dim} table. The available groups are {groups_avail}"
            )

    # Compute mean pw cc on all data layers - resulting in a df for each
    corr_dfs = [
        mean_pw_cc_by(ds, by, data_table=data_table, dim=dim)
        for data_table in data_tables
    ]

    # return a single merged df containing info for all data layer pw cc
    return reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        corr_dfs,
    )


def mean_pw_cc_by(ds, by, data_table="counts", dim="sample"):

    """Computes pairwise cc for all
    dim in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of dims
    in the group."""

    data = ds[f"{data_table}"].to_pandas()

    groups, pw_cc, n = [], [], []

    for s_group, group_ds in iter_groups(ds, by, dim):
        groups.append(s_group)
        n.append(len(group_ds[f"{dim}_id"].values))

        if len(group_ds[f"{dim}_id"].values) < 2:
            pw_cc.append(1.0)
            continue

        correlations = []
        for dim_ids in itertools.combinations(group_ds[f"{dim}_id"].values, 2):
            dim_0_enrichment = data.loc[:, dim_ids[0]]
            dim_1_enrichment = data.loc[:, dim_ids[1]]
            correlation = (
                st.pearsonr(dim_0_enrichment, dim_1_enrichment)[0]
                if np.any(dim_0_enrichment != dim_1_enrichment)
                else 1.0
            )
            correlations.append(correlation)
        pw_cc.append(round(sum(correlations) / len(correlations), 5))

    name = "_".join(by)
    column_prefix = f"{name}_{data_table}"

    ret = pd.DataFrame(
        {
            f"{name}": groups,
            f"{column_prefix}_pw_cc": np.array(pw_cc).astype(np.float64),
            f"{column_prefix}_n_reps": np.array(n).astype(np.int),
        }
    ).set_index(name)

    return ret


def collapse_sample_groups(*args, **kwargs):
    """wrap for sample collapse"""
    return collapse_groups(*args, **kwargs, collapse_dim="sample")


def collapse_peptide_groups(*args, **kwargs):
    """wrap for peptide collapse"""
    return collapse_groups(*args, **kwargs, collapse_dim="peptide")


def collapse_groups(
    ds, by, collapse_dim="sample", agg_func=np.mean, compute_pw_cc=False, **kwargs
):
    """
    Collapse a phip xarray dataset by some group in the metadata.
    """

    # Check to see if the group(s) is/are available
    groups_avail = ds[f"{collapse_dim}_metadata"].values
    for group in by:
        if group not in groups_avail:
            raise KeyError(
                f"{group} is not included as a column in the {collapse_dim} table. The available groups are {groups_avail}"
            )

    # define collapse and fixed df
    dims = set(["sample", "peptide"])
    fixed_dim = list(dims - set([collapse_dim]))[0]

    # grab relavent annotation tables
    collapse_df = ds[f"{collapse_dim}_table"].to_pandas()
    fixed_df = ds[f"{fixed_dim}_table"].to_pandas()

    # Create group-able dataset by assigning table columns to a coordinate
    if len(by) == 1:
        coord = collapse_df[by[0]]
        coord_ds = ds.assign_coords({f"{by[0]}": (f"{collapse_dim}_id", coord)})
    else:
        print(f"WARNING: Nothing available, here")
        return None

    # if were grouping by multiple things, we need to zip 'em into a tuple coord
    # psuedo-code
    # else:
    #    common_dim = f"{collapse_dim}_id"
    #    coor_arr = np.empty(len(ds[common_dim]), dtype=object)
    #    coor_arr[:] = list(zip(*(collapse_df[f].values for f in by)))
    #    coord_ds = ds.assign_coords(
    #        coord=xr.DataArray(coor_arr, collapse_dims=common_dim)
    #    )

    del coord_ds["sample_table"]
    del coord_ds["peptide_table"]
    del coord_ds["sample_metadata"]
    del coord_ds["peptide_metadata"]

    # Perform the reduction on all data tables.
    collapsed_enrichments = coord_ds.groupby(f"{by[0]}", squeeze=False).reduce(agg_func)

    if collapse_dim == "sample":
        collapsed_enrichments = collapsed_enrichments.transpose()

    # Once the data tables are grouped we have no use for first copy.
    del coord_ds

    # Compile each of the collapsed xarray variables.
    collapsed_xr_dfs = {
        f"{dt}": (
            ["peptide_id", "sample_id"],
            collapsed_enrichments[f"{dt}"].to_pandas(),
        )
        for dt in set(list(collapsed_enrichments.data_vars))
    }

    cat = throw_mismatched_features(collapse_df, by)

    # Compute mean pairwise correlation for all groups,
    # for all enrichment layers - and add it to the
    # resulting collapsed sample table
    if compute_pw_cc:
        mean_pw_cc = mean_pw_cc_by(ds, by, **kwargs)
        cat = cat.merge(mean_pw_cc, left_index=True, right_index=True)

    # Insert the correct annotation tables to out dictionary of variables
    collapsed_xr_dfs[f"{collapse_dim}_table"] = (
        [f"{collapse_dim}_id", f"{collapse_dim}_metadata"],
        cat,
    )

    collapsed_xr_dfs[f"{fixed_dim}_table"] = (
        [f"{fixed_dim}_id", f"{fixed_dim}_metadata"],
        fixed_df,
    )

    pds = xr.Dataset(
        collapsed_xr_dfs,
        coords={
            f"{collapse_dim}_id": collapsed_xr_dfs[f"{collapse_dim}_table"][
                1
            ].index.values,
            f"{fixed_dim}_id": collapsed_xr_dfs[f"{fixed_dim}_table"][1].index.values,
            f"{collapse_dim}_metadata": collapsed_xr_dfs[f"{collapse_dim}_table"][
                1
            ].columns.values,
            f"{fixed_dim}_metadata": collapsed_xr_dfs[f"{fixed_dim}_table"][
                1
            ].columns.values,
        },
    )
    return pds





