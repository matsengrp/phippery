r"""
=================
Escape Profile
=================

Use 
`optimal transport method <https://en.wikipedia.org/wiki/Transportation_theory_(mathematics)>`_
to compare 
`phage-dms <https://www.sciencedirect.com/science/article/pii/S2589004220308142>`_
escape profiles.
See the escape profile :ref:`description and examples <sec_escape_profile_comparisons>`.
"""

# dependencies
import pandas as pd
import numpy as np
import xarray as xr
import ot
from phippery.utils import peptide_id_coordinate_from_query
from phippery.utils import sample_id_coordinate_from_query
from Bio.Align import substitution_matrices

# TODO K: Start "helper functions" with _ i.e. the ones you dont want in documentation
# TODO k: doc strings of non "_" functions in RST format.

def get_aa_ordered_list():
    """
    return the ordered list of amino acid.
    This convention is based on the BLOSUM
    matrix in biopython and assumed for the
    binned distribution presenting amino
    acid contribution to differential 
    selection at a site
    """
    return list('ARNDCQEGHILKMFPSTWYV')


def get_cost_matrix():
    """
    return the default 40x40 cost matrix based on BLOSUM62
    and assigns maximum cost to transport between
    opposite signed differential selection contributions
    """

    substitution_matrix = substitution_matrices.load('BLOSUM62')
    alphabet_list = get_aa_ordered_list()
    Naa = len(alphabet_list)

    # chosen so that range of costs in the 
    # matrix is within an order of magnitude
    nthroot=7.

    # maximum cost assigned by the cost matrix
    maxMij = np.exp(np.max(-substitution_matrix)/nthroot)


    cost_matrix=[]

    # first 20 rows
    for aa in alphabet_list:
        row = [-x/nthroot for x in substitution_matrix[aa,:][:Naa]]
        cost_row = (np.exp(row)).tolist() + [maxMij for i in range(Naa)]
        cost_matrix.append(cost_row)

    # last 20 rows
    for aa in alphabet_list:
        row = [-x/nthroot for x in substitution_matrix[aa,:][:Naa]]
        cost_row = [maxMij for i in range(Naa)] + (np.exp(row)).tolist()
        cost_matrix.append(cost_row)

    return cost_matrix


def get_loc_esc_distr(
    ds,
    metric,
    sample_factor,
    sfact_val,
    loc
):
    """
    return the normalized distribution represented as a list
    for the amino acid pattern of scaled differential
    selection for a specified site and individual
    
    metric: label of the scaled differential selection data in ds
    loc: peptide annotation label for the location
    The individual is specified by a sample annotation label
    in sample_factor (e.g. 'sample_ID') and the corresponding
    value in sfact_val

    Parameters
    ----------

    ds : xarray.DataSet
        The dataset you would like to fit to
    """

    # Code assumes peptide annotation for location is called 'Loc',
    # and for amino acid substitution is called 'aa_sub'
    my_ds = ds.loc[
                dict(
                    peptide_id=peptide_id_coordinate_from_query(ds, [f"Loc == '{loc}'"]),
                    sample_id=sample_id_coordinate_from_query(ds, [f"{sample_factor} == '{sfact_val}'"]),
                )
            ]

    diff_sel = my_ds[metric].to_pandas().to_numpy().flatten()

    my_df = my_ds.peptide_table.loc[:,['aa_sub']].to_pandas()
    my_df['diff_sel'] = diff_sel

    esc_data_neg=[]
    esc_data_pos=[]
    for aa in get_aa_ordered_list():
        val = my_df[my_df['aa_sub']==aa]['diff_sel'].item()
        if val>0:
            esc_data_neg.append(0)
            esc_data_pos.append(val)
        else:
            esc_data_neg.append(-val)
            esc_data_pos.append(0)

    esc_data = esc_data_neg + esc_data_pos

    if np.sum(esc_data)==0:
        return esc_data
    else:
        return esc_data/np.sum(esc_data)


def get_weights(
    ds,
    metric,
    sample_factor,
    sfact_val1,
    sfact_val2,
    loc_start,
    loc_end
):
    """
    return the list of weights for computing the
    weighted sum of similarity scores in region
    [loc_start, loc_end]
    """

    # Code assumes peptide annotation for location is called 'Loc'

    loc_sums1=[]
    loc_sums2=[]
    for loc in range(loc_start, loc_end+1):
        ds1 = ds.loc[
                dict(
                    peptide_id=peptide_id_coordinate_from_query(ds, [f"Loc == '{loc}'"]),
                    sample_id=sample_id_coordinate_from_query(ds, [f"{sample_factor} == '{sfact_val1}'"]),
                    )
                ]

        diff_sel1 = ds1[metric].to_pandas().to_numpy().flatten()
        loc_sums1.append(0)
        for val in diff_sel1:
            loc_sums1[-1] = loc_sums1[-1] + abs(val)

        ds2 = ds.loc[
                dict(
                    peptide_id=peptide_id_coordinate_from_query(ds, [f"Loc == '{loc}'"]),
                    sample_id=sample_id_coordinate_from_query(ds, [f"{sample_factor} == '{sfact_val2}'"]),
                    )
                ]

        diff_sel2 = ds2[metric].to_pandas().to_numpy().flatten()
        loc_sums2.append(0)
        for val in diff_sel2:
            loc_sums2[-1] = loc_sums2[-1] + abs(val)

    loc_sums1 = loc_sums1/np.sum(loc_sums1)
    loc_sums2 = loc_sums2/np.sum(loc_sums2)

    weights={}
    total=0
    for i,loc in zip(range(loc_end-loc_start+1), range(loc_start, loc_end+1)):
        val = min(loc_sums1[i], loc_sums2[i])
        total = total+val
        weights[loc] = val

    weights = {k: v/total for k,v in weights.items()}

    return weights


def compute_sim_score(a, b, cost_matrix):
    """
    return the similarity score given
    two distributions and the cost matrix
    """

    if np.sum(a)==0 or np.sum(b)==0:
        return 0

    cost = ot.emd2(a, b, cost_matrix)
    return 1./cost


def loc_sim_score(
    ds,
    metric,
    cost_matrix,
    sample_factor,
    sfact_val1,
    sfact_val2,
    loc
):
    """
    return the similarity score for comparison at a site
    """

    a = get_loc_esc_distr(ds,metric,sample_factor,sfact_val1,loc)
    b = get_loc_esc_distr(ds,metric,sample_factor,sfact_val2,loc)

    return compute_sim_score(a, b, cost_matrix)


def region_sim_score(
    ds, 
    metric,
    cost_matrix,
    sample_factor,
    sfact_val1, 
    sfact_val2,
    loc_start,
    loc_end
):
    """
    return the similarity score for comparison in the 
    region [loc_start, loc_end]
    """

    weights = get_weights(ds, metric,
                          sample_factor, sfact_val1, sfact_val2,
                          loc_start, loc_end)
    region_sim=0
    for loc in range(loc_start, loc_end+1):
        sim = weights[loc] * loc_sim_score(ds, metric, cost_matrix,
                                           sample_factor, sfact_val1, sfact_val2,
                                           loc)
        region_sim = region_sim + sim

    return region_sim
