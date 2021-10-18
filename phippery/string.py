"""
@File: utils.py

@Author: Jared Galloway

Some helpful functions for printing
"""


# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import scipy.stats as st

# built-in python3
import os
import copy
import itertools
import io
from functools import reduce
from collections import defaultdict

from phippery.phipdata import get_sample_table
from phippery.phipdata import get_peptide_table
from phippery.phipdata import get_annotation_table


def string_feature(ds, feature: str, verbosity = 0, dim="sample"):
    """
    Take in a single feature and returna formatted string
    that gives a basic descriptiton of that feature and it's values.

    """
    t = get_sample_table(ds) if dim == "sample" else get_peptide_table(ds)
    ser = t[feature]
    dt = ser.dtype
    descript = """"""

    non_null_ser = ser[ser.notnull()]
    null_flag = len(non_null_ser) < len(ser)
    levels = list(set(non_null_ser.values))
    num_levels = len(levels)

    if dt == pd.StringDtype():

        if num_levels >= 50:
            descript += f"""
Sorry, {feature} is a String feature with too many factor levels 
(n = {num_levels}) to describe, here 
"""

        elif num_levels < 50 and num_levels >=2:
            descript += f"""
{feature}: String Feature
-------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

> {feature} in ['{levels[0]}', '{levels[1]}', ...]

> {feature} not in ['{levels[0]}', '{levels[-1]}', ...]

> {feature} != '{levels[-2]}'
"""

        else:
            descript += f"""
There's only a single factor level, {levels[0]}, across all samples.
"""

    elif dt == pd.BooleanDtype():
        descript += f"""
{feature}: Boolean Feature:
---------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

> {feature} == True

> {feature} == False
"""

    elif dt == pd.Int64Dtype():
        des = ser.describe()
        descript += f"""
{feature}: Integer Feature:
---------------------------

distribution of numerical feature:

{des}


Some example query statements:
------------------------------

> {feature} >= {int(des[1])}

> {feature} <= {int(des[1])}

> ({feature} >= {int(des[4])}) and ({feature} <= {int(des[5])})
"""

    elif  dt == pd.Float64Dtype():

        des = ser.describe()
        descript += f"""
{feature}: Float Feature:
-------------------------

distribution of numerical feature:

{des}


Some example query statements:
------------------------------

> {feature} >= {des[1]}

> {feature} <= {des[1]}

> ({feature} >= {des[4]}) and ({feature} <= {des[5]})
"""

    else:
        st.error("Never seen this Dtype before!")

    if null_flag:
        descript += f"""
> {feature}.isnull()

> {feature}.notnull()
"""

    return descript


def string_ds(ds, verbosity: int):
    """
    Summarize the data in a given dataset

    If verbosity flag is set to zero, this will print the 
    basic information about number of
    samples, peptides, and enrichment layers
    in a given dataset. With a verbosity of one (-v) you will
    get also information about annotations and available datatypes.
    If verbosity flag is set to two (-vv) - Print 
    detailed information about all data tables including annotation
    feature types and distribution of enrichment values for each layer.
    A verbosity of three will basically loop through all features
    """

    # initialize formatting for each of the three major facets of a dataset
    table_strings = {
        'sample_table' : """
Sample Table:
-------------\n""",

        'peptide_table' : """
Peptide Table:
--------------\n""",

        'enrichments' : """
Enrichment Matrices:
--------------------\n"""
    }

    for dimension in ["sample", "peptide"]:

        #num_dimensions= len(ds[f"{dimension}_id"].values)
        #dimension_annotations = list(ds[f"{dimension}_metadata"].values)

        
        #dimension_annotation_strings = {
        #    f'{l}':f"""\n           * {l}"""
        #for l in dimension_annotations
        #}

        ## call
        #if verbosity > 0:
        #    pass
        #if verbosity > 1:
        #    pass

        #complete = """"""
        #for key, value in dimension_annotation_strings.items():
        #    complete += value

        df = get_annotation_table(ds, dim=dimension)
        num_dimensions = len(df) 

        buffer = io.StringIO()
        df.info(buf=buffer)
        complete = buffer.getvalue()
        #* Number of {dimension}s: {num_dimensions}

        #* {dimension} annotation table features: 
        
        table_strings[f"{dimension}_table"] += f"""
        {complete}
        """

    # initialize formatting strings for all enrichment layers
    enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    enrichment_strings = {}
    for enr in enr_layers:
        mat = ds[enr].to_pandas()
        enrichment_strings[enr] = f"""* {enr}\n{mat.describe()}"""

    #enrichment_strings = {
    #    f'{l}':f"""             * {l}"""
    #for l in enr_layers
    #}

    complete = """"""
    for key, value in enrichment_strings.items():
        complete += value
    table_strings['enrichments'] += complete

    final = """"""
    for key, value in table_strings.items():
        final += value

    return final
