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

from phippery.utils import get_annotation_table


def string_feature(ds, feature: str, verbosity=0, dim="sample", numeric_dis=True):
    """
    Take in a single feature and returna formatted string
    that gives a basic descriptiton of that feature and it's values.

    """
    t = get_annotation_table(ds, dim)
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

        elif num_levels < 50 and num_levels >= 2:
            descript += f"""
{feature}: {dt} Feature
-------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

{feature} in ['{levels[0]}', '{levels[1]}', ...]

{feature} not in ['{levels[0]}', '{levels[-1]}', ...]

{feature} != '{levels[-2]}'
"""

        else:
            descript += f"""
There's only a single factor level, {levels[0]}, across all samples.
"""

    elif dt == pd.BooleanDtype():
        descript += f"""
{feature}: {dt} Feature:
---------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

{feature} == True

{feature} == False
"""

    elif dt == pd.Int64Dtype() or dt == pd.Float64Dtype():
        if num_levels == 1:
            descript += f"""
There's only a single factor level, {levels[0]}, across all samples.
"""

        elif (num_levels > 1) and (not numeric_dis):
            descript += f"""
{feature}: {dt} Feature
-------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

{feature} in [{levels[0]}, {levels[1]}, ...]

{feature} not in [{levels[0]}, {levels[-1]}, ...]

{feature} != {levels[-2]}
"""
        else:
            des = ser.describe()
            descript += f"""
{feature}: {dt} Feature:
---------------------------

distribution of numerical feature:

{des}

Some example query statements:
------------------------------

{feature} >= {int(des[1])}

{feature} <= {int(des[1])}

({feature} >= {int(des[4])}) and ({feature} <= {int(des[5])})
"""

    else:
        descript += f"""
{feature}: {dt} Feature
-------------------------

Factor level counts df:

{non_null_ser.value_counts()}


Some example query statements:
------------------------------

{feature} in ['{levels[0]}', '{levels[1]}', ...]

{feature} not in ['{levels[0]}', '{levels[-1]}', ...]

{feature} != '{levels[-2]}'
"""

    if null_flag:
        descript += f"""
{feature}.isnull()

{feature}.notnull()
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
        "sample_table": """
Sample Table:
-------------\n""",
        "peptide_table": """
Peptide Table:
--------------\n""",
        #        'enrichments' : """
        # Enrichment Matrices:
        # --------------------\n"""
    }

    for dimension in ["sample", "peptide"]:

        df = get_annotation_table(ds, dim=dimension)
        num_dimensions = len(df)

        buffer = io.StringIO()
        df.info(buf=buffer)
        complete = buffer.getvalue()
        table_strings[f"{dimension}_table"] += f"""{complete}"""

    # initialize formatting strings for all enrichment layers
    # enr_layers = set(list(ds.data_vars)) - set(["sample_table", "peptide_table"])
    # enrichment_strings = {}
    # for enr in enr_layers:
    #    mat = ds[enr].to_pandas()
    #    enrichment_strings[enr] = f"""* {enr}\n{mat.describe()}"""

    # complete = """"""
    # for key, value in enrichment_strings.items():
    #    complete += value
    # table_strings['enrichments'] += complete

    final = """"""
    for key, value in table_strings.items():
        final += value

    return final
