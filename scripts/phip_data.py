"""
file: phip_data.py

@author: Jared Galloway

Some simple functions to get us thinking about 
useful operations we can do to the data.
"""
import pandas as pd
import argparse
import scipy.stats as st
import re
import sys
import os
import copy
from matplotlib import pyplot as plt

def extract_sample_info(
    raw_counts):

    # TODO there's certainly a cleaner way to do this.
    technical_rep_1 = [re.match("\d+\.[1]",f"{sample_rep}")[0] 
        for sample_rep in raw_counts.columns 
        if re.match("\d+\.[1]",f"{sample_rep}") != None]

    technical_rep_2 = [re.match("\d+\.[2]",f"{sample_rep}")[0] 
        for sample_rep in raw_counts.columns 
        if re.match("\d+\.[2]",f"{sample_rep}") != None]

    return technical_rep_1, technical_rep_2


def prune_non_correlated_replicate_samples(
    raw_counts, 
    technical_rep_1,
    technical_rep_2,
    corr_thresh = 0.9, 
    pval_thresh = 0.05,
    beads_only = '35'):

    pruned_counts = raw_counts.copy()

    correlations, samples = [], []
    for i, (t1, t2) in enumerate(zip(technical_rep_1, technical_rep_2)):
        
        # let's make sure we're looking at the same sample
        assert(t1.split(".")[0] == t2.split(".")[0])
        sample = t1.split(".")[0]

        # Now, lets determine which samples to drop,
        # and which samples we can sum replicates over
        # based off of correlation.
        corr = st.pearsonr(pruned_counts[t1], pruned_counts[t2])
        correlations.append(corr)
        samples.append(sample)
            
        # if the technical replicates are correlated enough then
        # sum the pruned_counts, add 10, and insert new column.
        #print(f"sample {sample} has technical rep correlation: {corr[0]} \
        #    with a p-value of {round(corr[1],5)}")
        if (corr[0] > corr_thresh and corr[1] <= pval_thresh) or sample == beads_only:
            pruned_counts[sample] = pruned_counts[t1] + pruned_counts[t2] + 20
        #else:
        #    print(f"dropping sample {sample}, b/c corr was: {corr[0]}")

        pruned_counts.drop([t1,t2], axis=1, inplace=True)
        
    return pruned_counts, correlations, samples



def get_standardized_enrichment(
    raw_counts, 
    beads_only_sample = "35", 
    library_control_sample = "37.input"):

    """
    Standardize a counts matrix by 
        1. converting each column (sample) into frequencies
        2. dividing each item in row by the total abundance of
            phage ina library (37.input)
        3. subtract background (negative control / mock IP) 
            from background.

    A small example is below.

    >>> import numpy as np
    >>> import pandas as pd
    >>> np.random.seed(23)
    >>> 
    >>> df = pd.DataFrame(np.random.randint(10, size=12).reshape(3,4),columns=list('abcd'))
    >>> df
       a  b  c  d
    0  3  6  8  9
    1  6  8  7  9
    2  3  6  1  2
    >>> freq = df.div(df.sum(axis=0), axis=1)
    >>> freq
          a    b       c     d
    0  0.25  0.3  0.5000  0.45
    1  0.50  0.4  0.4375  0.45
    2  0.25  0.3  0.0625  0.10
    >>> enr = freq.div(freq["d"], axis=0)
    >>> enr
              a         b         c    d
    0  0.555556  0.666667  1.111111  1.0
    1  1.111111  0.888889  0.972222  1.0
    2  2.500000  3.000000  0.625000  1.0
    >>> enr["e"] = [1, 2, 3]
    >>> enr
              a         b         c    d  e
    0  0.555556  0.666667  1.111111  1.0  1
    1  1.111111  0.888889  0.972222  1.0  2
    2  2.500000  3.000000  0.625000  1.0  3
    >>> std_enr = enr.sub
    enr.sub(       enr.subtract(  
    >>> std_enr = enr.subtract(enr["e"], axis=0)
    >>> std_enr
              a             b         c    d  e
    0 -0.444444 -3.333333e-01  0.111111  0.0  0
    1 -0.888889 -1.111111e+00 -1.027778 -1.0  0
    2 -0.500000 -4.440892e-16 -2.375000 -2.0  0
    """
    frequencies = raw_counts.truediv(raw_counts.sum(axis=0), axis=1)
    enrichment = frequencies.truediv(frequencies[library_control_sample], axis=0)
    standardized_enrichment = enrichment.subtract(enrichment[beads_only_sample], axis=0)

    return standardized_enrichment






