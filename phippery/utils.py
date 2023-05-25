r"""
=================
Utils
=================

Utilities for building, indexing, and manipulating
and xarray dataset topology 
specific to most **phippery** functions provided in this package
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


def iter_groups(ds, by, dim="sample"):
    """This function returns an iterator
    yeilding subsets of the provided dataset,
    grouped by items in the metadata of either of the
    dimensions specified.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset to iterate over.

    Returns
    -------

    generator :
        Returns subsets of the original dataset sliced by
        either sample or peptide table groups.

    Example
    -------

    >>> phippery.get_annotation_table(ds, "sample")
    sample_metadata  fastq_filename reference seq_dir sample_type
    sample_id
    0                sample_0.fastq      refa    expa  beads_only
    1                sample_1.fastq      refa    expa  beads_only
    2                sample_2.fastq      refa    expa     library
    3                sample_3.fastq      refa    expa     library    
    >>> ds["counts"].values
    array([[458, 204, 897, 419],
           [599, 292, 436, 186],
           [ 75,  90, 978, 471],
           [872,  33, 108, 505],
           [206, 107, 981, 208]])
    >>> sample_groups = iter_groups(ds, by="sample_type")
    >>> for group, phip_dataset in sample_groups:
    ...     print(group)
    ...     print(phip_dataset["counts"].values)
    ...
    beads_only
    [[458 204]
     [599 292]
     [ 75  90]
     [872  33]
     [206 107]]
    library
    [[897 419]
     [436 186]
     [978 471]
     [108 505]
     [981 208]]
    """

    table = get_annotation_table(ds, dim=dim)
    for group, group_df in table.groupby(by):
        group_ds = ds.loc[{f"{dim}_id": list(group_df.index.values)}]
        yield group, group_ds


def get_annotation_table(ds, dim="sample"):
    """
    return a copy of the peptide table after converting all
    the data types applying the pandas NaN heuristic

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset to extract an annotation from.

    dim : str
        The annotation table to grab: "sample" or "peptide".

    Returns
    -------

    pd.DataFrame :
        The annotation table.

    """

    st = ds[f"{dim}_table"].to_pandas().convert_dtypes()
    st.index.name = f"{dim}_id"
    return st


def stitch_dataset(
    counts,
    peptide_table,
    sample_table,
):
    """Build an phippery xarray dataset from individual
    tables.

    Note
    ----
    If the counts matrix that you're passing has the shape
    (M x N) for M peptides, and N samples, the
    sample table should have a len of N,
    and peptide table should have len M

    Parameters
    ----------

    counts : numpy.ndarray
        The counts matrix for sample peptide enrichments.

    sample_table : pd.DataFrame
        The sample annotations corresponding to the columns of
        the counts matrix.

    peptide_table : pd.DataFrame
        The peptide annotations corresponding to the rows of
        the counts matrix.

    Returns
    -------

    xarray.DataSet :
        The formatted phippery xarray dataset used by most of the phippery functionality.

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


def collect_counts(counts):
    r"""merge individual tsv files for a bunh of samples
    into a counts matrix.

    Parameters
    ----------

    counts : list[str]
        A list a filepaths relative to current working directory to read in.
        The filepaths should point to tab-separated files for each sample
        which contains two columns (without headers):

            1. peptide ids - the integer peptide identifiers
            2. enrichments - the respective enrichments for any peptide id

    Returns
    -------

    pd.DataFrame :
        The merged enrichments with peptides as the index, and filenames as column names.

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


def to_tall(ds: xr.Dataset):
    """Melt a phippery xarray dataset into a single long-formatted
    dataframe that has a unique sample peptide interaction on each
    row. Ideal for ggplotting.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset to extract an annotation from.

    Returns
    -------

    pd.DataFrame :
        The tall formatted dataset.

    Example
    -------
    >>> ds["counts"].to_pandas()
    sample_id     0    1
    peptide_id
    0           453  393
    1           456  532
    2           609  145
    >>> to_tall(ds)[["sample_id", "peptide_id", "counts"]]
      sample_id  peptide_id  counts
    0         0           0     453
    1         0           1     456
    2         0           2     609
    3         1           0     393
    4         1           1     532
    5         1           2     145
    """

    return pd.concat([
        sample_df
        for sample_df in yield_tall(ds)
    ])

def yield_tall(ds: xr.Dataset):
    """For each sample, yield a tall DataFrame."""

    # Get the table of samples
    sample_table = ds.sample_table.to_pandas().reset_index().infer_objects()

    # Keep track of the order of columns in all emitted items
    cnames = None

    # For each sample
    for sample_id in sample_table["sample_id"].values:

        # Make sure that values for this sample are present in all data tables
        for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"]):
            assert sample_id in ds[f"{dt}"].to_pandas().columns.values, f"Could not find sample '{sample_id}' in table for {dt}"

        # Make a wide table
        sample_df = pd.DataFrame({
            f"{dt}": ds[
                f"{dt}"
            ].to_pandas(
            )[
                sample_id
            ]
            for dt in set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
        }).assign(
            sample_id=sample_id
        )

        # Get the table of peptides
        peptide_table = ds.peptide_table.to_pandas().reset_index().infer_objects()

        # merge the metadata into the melted datatables
        sample_df = sample_df.merge(peptide_table, on="peptide_id")
        
        # Merge the sample table
        sample_df = sample_df.merge(sample_table, on="sample_id")

        # If the column names have not yet been set
        if cnames is None:

            # Set them based on the first table
            cnames = sample_df.columns.values

        # If the column names were set in a previous iteration
        else:

            # Make sure that this table conforms to the same column order
            sample_df = sample_df.reindex(columns=cnames)

        yield sample_df


def to_wide(ds):
    """Take a phippery dataset and split it into
    its separate components in a dictionary.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset to separate.

    Returns
    -------
    dict :
        The dictionary of annotation tables and enrichment matrices.

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


def to_wide_csv(ds, file_prefix):
    """Take a phippery dataset and split it into
    its separate components at writes each into a
    comma separated file.

    Note
    ----
    This is the inverse operation of the
    `dataset_from_csv()` utility function.
    Generally speaking these functions are used for 
    long term storage in common formats when pickle
    dumped binaries are not ideal.


    Parameters
    ----------

    ds : xarray.DataSet
        The dataset to extract an annotation from.

    file_prefix : str
        The fileprefix relative to the current working directory
        where the files should be written.

    Returns
    -------
    None
    """

    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    for dt in enr_layers:
        layer = ds[f"{dt}"].to_pandas()
        layer.index.name = ""
        layer.to_csv(f"{file_prefix}_{dt}.csv", na_rep="NA", index_label=None)

    for at in ["sample", "peptide"]:
        get_annotation_table(ds, dim=at).to_csv(
            f"{file_prefix}_{at}_annotation_table.csv"
        )


def id_coordinate_from_query_df(ds, query_df):
    """Given a dataframe with pandas query statements
    for both samples and peptides, return the relevent sample
    and peptide id's after applying the logical AND of all queries.

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to query.

    query_df : pd.DataFrame
        A dataframe with must have two columns (including headers):

            1. "dimension" - either "sample" or "peptide" to specify expression dimension
            2. "expression" - The pandas query expression to apply.

    Returns
    -------
        tuple : list, list
            Return a tuple of sample id's and peptide id's

    """

    sq = list(query_df.loc[query_df["dimension"] == "sample", "expression"].values)
    sid = id_query(ds, "sample", " & ".join(sq))

    pq = list(query_df.loc[query_df["dimension"] == "peptide", "expression"].values)
    pid = id_query(ds, "peptide", " & ".join(pq))

    return sid, pid


def ds_query(ds, query, dim="sample"):
    """
    Apply a sample or peptide query statement to the
    entire dataset.

    Note
    ----
    For more on pandas queries, see
    `the pandas documentation <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to query.

    query : str
        pandas query expression

    dim : str
        The dimension to to apply the expression

    Returns
    -------
    xarray.DataSet :
        reference to the dataset slice from the given expression.

    Example
    -------

    >>> phippery.get_annotation_table(ds, "peptide")
    peptide_metadata Oligo   virus
    peptide_id
    0                 ATCG    zika
    1                 ATCG    zika
    2                 ATCG    zika
    3                 ATCG    zika
    4                 ATCG  dengue
    5                 ATCG  dengue
    6                 ATCG  dengue
    7                 ATCG  dengue
    >>> zka_ds = ds_query(ds, "virus == 'zika'", dim="peptide")                                              
    >>> zka_ds["counts"].to_pandas()                                                                         
    sample_id     0    1    2    3    4    5    6    7    8    9
    peptide_id
    0           110  829  872  475  716  815  308  647  216  791
    1           604  987  776  923  858  985  396  539   32  600
    2           865  161  413  760  422  297  639  786  857  878
    3           992  354  825  535  440  416  572  988  763  841
    """

    idx = id_query(ds, query, dim)
    return ds.loc[{f"{dim}_id": idx}]


def id_query(ds, query, dim="sample"):
    """
    Apply a sample or peptide query statement to the
    entire dataset and retrieve the respective indices.

    Note
    ----
    For more on pandas queries, see
    https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to query.

    query : str
        pandas query expression

    dim : str


    Return
    ------
    list[int] :
        The list of integer identifiers that apply to a given expression
        for the respective dimension.
    """

    return get_annotation_table(ds, dim).query(query).index.values


def load(path):
    """simple wrapper for loading
    xarray datasets from pickle binary

    Parameters
    ----------

    path : str
        Relative path of binary phippery dataset

    Returns
    -------

    xarray.DataSet :
        phippery dataset
    """

    return pickle.load(open(path, "rb"))


def dump(ds, path):
    """
    simple wrapper for dumping xarray datasets to pickle binary

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to write to disk.

    path : str
        The relative path you would like to write to.

    Returns
    -------
    None
    """

    pickle.dump(ds, open(path, "wb"))


def dataset_from_csv(
    peptide_table_filename, 
    sample_table_filename, 
    counts_table_filename
):
    r"""Load a dataset from individual comma separated
    files containing the counts matrix, as well as
    sample and peptide annotation tables.

    Note
    ----
    This is the inverse operation of the
    `to_wide_csv()` utility function.
    Generally speaking these functions are used for 
    long term storage in common formats when pickle
    dumped binaries are not ideal.
    For now, this function only supports
    a single enrichment table to be added with the
    variable name "counts" to the dataset.
    If you would like to add other transformation of the
    enrichment table (i.e. cpm, mlxp, etc), you can
    load the csv's via pandas and add to the dataset
    using the `add_enrichment_layer_from_array` function


    Parameters
    ----------
    counts_table_filename : str
        The glob filepath to csv file(s) containing the enrichments.
        All files should have the first column be indices which match the
        given peptide table index column.
        The first row then should have column headers that match the
        index of the sample table.

    peptide_table_filename : str
        The relative filepath to the peptide annotation table.

    sample_table_filename : str
        The relative filepath to the sample annotation table.

    Returns
    -------
    xarray.DataSet :
        The combined tables in a phippery dataset.
    """

    counts_df = _collect_counts_matrix(counts_table_filename)
    peptide_df = _collect_peptide_table(peptide_table_filename)
    sample_df = _collect_sample_table(sample_table_filename)

    return stitch_dataset(
        counts=counts_df,
        peptide_table=peptide_df,
        sample_table=sample_df,
    )


def _collect_sample_table(sample_table_filename: str):
    """Read and verify a sample table."""

    sample_table = pd.read_csv(sample_table_filename, sep=",", index_col=0, header=0)

    if sample_table.index.name != "sample_id":
        raise ValueError("The name of the index must be 'sample_id'")

    if sample_table.index.dtype != "int64":
        raise ValueError("The index values for sample_id must be inferred as integers")

    sample_table.sort_index(inplace=True)
    return sample_table


def _collect_peptide_table(peptide_table_filename: str):
    """Read and verify a peptide table."""

    peptide_table = pd.read_csv(peptide_table_filename, sep=",", index_col=0, header=0)

    if peptide_table.index.name != "peptide_id":
        raise ValueError

    if peptide_table.index.dtype != "int64":
        raise ValueError("The index values for peptide_id must be inferred as integers")

    peptide_table.sort_index(inplace=True)
    return peptide_table


def _collect_counts_matrix(counts_matrix_filename: str):
    """Read and verify a counts matrix file."""

    counts_matrix = pd.read_csv(counts_matrix_filename, sep=",", index_col=0, header=0)

    try:
        counts_matrix.columns = counts_matrix.columns.astype(int)
    except:
        raise ValueError(
            "column header values much be able to cast to type 'int' to match peptide table index"
        )

    try:
        counts_matrix.index = counts_matrix.index.astype(int)
    except:
        raise ValueError(
            "row index values much be able to cast to type 'int' to match peptide table index"
        )

    counts_matrix.sort_index(inplace=True)
    counts_matrix.sort_index(axis=1, inplace=True)

    return counts_matrix


def add_enrichment_layer_from_array(ds, enrichment, new_table_name=None, inplace=True):
    """Append an enrichment layer to the dataset.

    Parameters
    ----------

    ds : xarray.DataSet
        The phippery dataset to append to.

    enrichment : np.array
        The enrichment matrix to append to the phippery dataset.
        The number of rows should be the same length as ds.peptide_id
        and the number of columns should be the same length as ds.sample_id

    new_table_name : str
        What you would like to name the enrichment layer.

    inplace : bool
        Determines whether to modify the passed dataset, or return an augmented
        copy

    Returns
    -------
    None | xarray.DataSet
        The augmented phippery dataset copy is returned if inplace is ``True``
    """

    if enrichment.shape != ds.counts.shape:
        ins = enrichment.shape
        cur = ds.counts.shape
        pri = f"provided enrichment layer shape: {ins},"
        pri += f"current working dataset counts shape: {cur}"
        raise ValueError(
            f"Enrichments must have the same shape as enrichments in dataset. {pri}"
        )
    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    if new_table_name == None:
        new_table_name = f"enrichment_layer_{len(enr_layers)+1}"
    if inplace:
        ds[new_table_name] = xr.DataArray(enrichment, dims=ds.counts.dims)
        return None
    else:
        ds_copy = copy.deepcopy(ds)
        ds_copy[new_table_name] = xr.DataArray(enrichment, dims=ds.counts.dims)
        return ds_copy


def _throw_mismatched_features(df, by):
    """
    When you collapse by some set of columns in the dataframe,
    keep only features which homogeneous within groups.

    This is similar to 'DataFrameGroupby.first()' method,
    but instead of keeping the first factor level appearing for each group
    category, we only throw any features which are note homogeneous within groups.
    You could achieve the same functionality by listing features you know to be
    homogeneous in the 'by' parameter.
    """

    # Find out which collapse features are shared within groups
    collapsed_sample_metadata = defaultdict(list)
    for i, (group, group_df) in enumerate(df.groupby(by)):
        for column, value in group_df.items():
            v = value.values
            if np.all(v == v[0]) or np.all([n != n for n in v]):
                collapsed_sample_metadata[column].append(v[0])

    # Throw out features that are not shared between groups
    to_throw = [
        key for key, value in collapsed_sample_metadata.items() if len(value) < i + 1
    ]
    [collapsed_sample_metadata.pop(key) for key in to_throw]
    return pd.DataFrame(collapsed_sample_metadata)


def _mean_pw_cc_by_multiple_tables(ds, by, dim="sample", data_tables="all"):
    """
    A wrapper for computing pw cc within groups defined
    with the 'by' parameter. Merges the correlations into a
    single table
    """

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
        _mean_pw_cc_by(ds, by, data_table=data_table, dim=dim)
        for data_table in data_tables
    ]

    # return a single merged df containing info for all data layer pw cc
    return reduce(
        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
        corr_dfs,
    )


def _mean_pw_cc_by(ds, by, data_table="counts", dim="sample"):
    """
    Computes pairwise cc for all
    dim in a group specified by 'group' column.

    returns a dataframe with each group, it's
    repective pw_cc, and the number of dims
    in the group.
    """

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
            f"{column_prefix}_n_reps": np.array(n).astype(int),
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
    """Collapse an xarray dataset along one of the annotation axis
    by applying the agg_function to annotation groups of 'by'.

    Parameters
    ----------

    ds : xarray.DataSet
        The phippery dataset to append to.

    by : list
        The name of the annotation feature you would like to collapse.

    collapse_dim : str
        The dimension you's like to collapse. "sample" or "peptide"

    compute_pw_cc : bool
        Whether or not to compute the mean pairwise correlation
        of all values within any feature group that is being
        collapsed.

    agg_func : *function*
        This function must take a one-dimensional array and aggregate
        all values to a single number, agg_func(list[float | int]) -> float | int


    Returns
    -------

    xarray.DataSet :
        The collapsed phippery dataset.

    Example
    -------
    >>> get_annotation_table(ds, dim="sample")
    sample_metadata  fastq_filename reference seq_dir sample_type
    sample_id
    0                sample_0.fastq      refa    expa  beads_only
    1                sample_1.fastq      refa    expa  beads_only
    2                sample_2.fastq      refa    expa     library
    3                sample_3.fastq      refa    expa     library
    4                sample_4.fastq      refa    expa          IP
    5                sample_5.fastq      refa    expa          IP
    >>> ds["counts"].to_pandas()
    sample_id   0  1  2  3  4  5
    peptide_id
    0           7  0  3  2  3  2
    1           6  3  1  0  7  5
    2           9  1  7  8  4  7
    >>> mean_sample_type_ds = collapse_groups(ds, by=["sample_type"])
    >>> get_annotation_table(mean_sample_type_ds, dim="sample")
    sample_metadata reference seq_dir sample_type
    sample_id
    0                    refa    expa          IP
    1                    refa    expa  beads_only
    2                    refa    expa     library
    >>> mean_sample_type_ds["counts"].to_pandas()
    sample_id     0    1    2
    peptide_id
    0           2.5  3.5  2.5
    1           6.0  4.5  0.5
    2           5.5  5.0  7.5
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

    cat = _throw_mismatched_features(collapse_df, by)

    # Compute mean pairwise correlation for all groups,
    # for all enrichment layers - and add it to the
    # resulting collapsed sample table
    if compute_pw_cc:
        mean_pw_cc = _mean_pw_cc_by(ds, by, **kwargs)
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
