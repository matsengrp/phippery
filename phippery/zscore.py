"""
@File: zscore.py

@Author: Kevin Sung
"""

import numpy as np
import pandas as pd
from functools import reduce

# Binning is represented as an array of arrays: an element array stores the 'peptide_id's in a bin
def zscore_pids_binning(
    beads_ds,  # dataset of beads-only samples
    data_table,  # peptide quantity for performing binning
    min_Npeptides_per_bin,  # mininum number of 'peptide_id's in a bin
):
    beads_data = beads_ds[data_table].values

    # For each peptide, sum the 'data_table' quantity across beads-only samples
    # Initially, all 'peptide_id's with the same sum are binned together
    sums_df = pd.DataFrame(index=beads_ds.peptide_id.values, columns=["sum"])
    sums_df["sum"] = np.sum(beads_data, axis=1)
    sums_df = sums_df.sort_values("sum")
    uniq_sum_values = sums_df["sum"].drop_duplicates().values

    binning = []
    for val in uniq_sum_values:
        val_pids = sums_df[sums_df["sum"] == val].index.tolist()
        binning.append(val_pids)

    # For bins with less than 'min_Npeptides_per_bin' peptides,
    # merge with adjacent bin of lower number of peptides
    done = False
    iprev = 0
    while done == False:
        done = True
        istart = iprev
        for i in range(istart, len(binning)):
            iprev = i
            if len(binning[i]) < min_Npeptides_per_bin:
                done = False

                # first bin in binning
                if i == 0:
                    if len(binning) == 1:  # if this is the only bin, we are done
                        done = True
                        break
                    else:
                        binning[i : i + 2] = [
                            reduce(lambda j, k: j + k, binning[i : i + 2])
                        ]

                # last bin in binning
                elif i == len(binning) - 1:
                    binning[i - 1 : i + 1] = [
                        reduce(lambda j, k: j + k, binning[i - 1 : i + 1])
                    ]

                else:  # merge with adjacent bin of lower number of peptides
                    if len(binning[i - 1]) > len(binning[i + 1]):
                        binning[i : i + 2] = [
                            reduce(lambda j, k: j + k, binning[i : i + 2])
                        ]
                    else:
                        binning[i - 1 : i + 1] = [
                            reduce(lambda j, k: j + k, binning[i - 1 : i + 1])
                        ]
                break

    return binning


# Compute Z-scores for each sample in the dataset
def compute_zscore(
    ds,
    data_table,  # peptide quantity for computing z-scores
    binning,  # binning of 'peptide_id's
    lower_quantile_limit,  # counts below this quantile are ignored for computing mean, stddev
    upper_quantile_limit,  # counts above this quantile are igonred for computing mean, stddev
):
    data_df = ds[data_table].to_pandas()

    # construct tables for means and stddevs in each bin (row)
    # across all samples in the dataset (column)
    means_df = pd.DataFrame(columns=data_df.columns)
    stds_df = pd.DataFrame(columns=data_df.columns)
    for pids in binning:
        bin_df = ds.loc[dict(peptide_id=pids)][data_table].to_pandas()
        bin_means = []
        bin_stds = []
        for col in bin_df.columns:
            bin_data = bin_df[col].to_numpy()
            lower_limit = np.quantile(bin_data, lower_quantile_limit)
            upper_limit = np.quantile(bin_data, upper_quantile_limit)
            null_data = bin_data[(bin_data >= lower_limit) & (bin_data <= upper_limit)]
            null_mean = np.mean(null_data)
            null_std = np.std(null_data)
            bin_means.append(null_mean)
            bin_stds.append(null_std)
        means_df.loc[len(means_df.index)] = bin_means
        stds_df.loc[len(stds_df.index)] = bin_stds

    # construct mapping of peptide_id to bin number
    pid_to_bin_df = pd.DataFrame(index=ds.peptide_id)
    pid_to_bin_df["ibin"] = -1
    for i in range(len(binning)):
        for pid in binning[i]:
            pid_to_bin_df.loc[pid] = i

    zscore_df = pd.DataFrame(index=data_df.index, columns=data_df.columns)
    for pid in zscore_df.index:
        ibin = pid_to_bin_df.loc[pid]["ibin"]
        if ibin < 0:
            print("peptide_id", pid, "is not in binning!")
            zs = 0
        else:
            data = data_df.loc[pid].to_numpy()
            null_means = means_df.loc[ibin].to_numpy()
            null_stds = stds_df.loc[ibin].to_numpy()
            zs = (data - null_means) / null_stds
        zscore_df.loc[pid] = zs

    # Convert Infs and NaNs to zero
    zscore_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    zscore_df.fillna(value=0, inplace=True)

    return zscore_df, means_df, stds_df
