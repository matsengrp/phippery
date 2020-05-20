"""
@File: fold_analysis

@Author: Jared Galloway

UNDER CONSTRUCTION

This includes functions that
compute fole enrichment analysis


* compile tsv files into a phip dataset
* TODO check phip_dataset attributed

"""

# dependencies
import pandas as pd
import numpy as np
import scipy.stats as st


def calculate_fold_normalized_enrichment(
    phip_dataset, library_control_sample, mock_ip_sample
):
    """
    calculate the fold enrichment counts matrix,
    currently, we are expecting a dictionary-type dataset
    for counts, peptide metadata, and sample metadata
    """
    frequencies = phip_dataset["counts"].truediv(
        phip_dataset["counts"].sum(axis=0), axis=1
    )
    # no inplace :(
    # TODO, think about how to save some memory here
    del phip_dataset["counts"]
    enrichment = frequencies.truediv(frequencies[library_control_sample], axis=0)
    del frequencies
    standardized_enrichment = enrichment.subtract(enrichment[mock_ip_sample], axis=0)
    phip_dataset["counts"] = standardized_enrichment

    return None
