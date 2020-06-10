"""
@File: plotting.py

@Author: Jared Galloway

This should include some helpful function which
will allow us to plot interesting things given
a PhipData object
"""

import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import cm
import re
from phippery.PhipData import *


def plot_peptide_enrichment_by_nt_position(
    ds, strain_pattern, sample, out, cmap="inferno"
):
    """
    This function takes plots a subset of tile enrichments
    as a function of nucleotide position (as described by
    nt_start, nt_end in peptide metadata) for all specific sample.
    # TODO bad practice to have a variable take on multiple types.
    :param: strain_pattern -
    if this parameter is a string, it is assumed to be a regex pattern
    and all virus strains in the ds that match are plotted. Otherwise,
    a list of viruses should be provided to
    """
    # TODO assert phip dataset consistency

    # if type(strain_pattern) == str:
    all_strains = set(ds["peptide_metadata"]["Virus_Strain"])
    selected_strains = [
        re.match(strain_pattern, f"{strain}")[0]
        for strain in all_strains
        if re.match(strain_pattern, f"{strain}") is not None
    ]
    assert len(all_strains) != 0

    fig, ax = plt.subplots()
    color_map = cm.get_cmap(cmap, len(selected_strains)).colors
    for i, strain in enumerate(selected_strains):

        strain_indices = ds["peptide_metadata"][
            ds["peptide_metadata"]["Virus_Strain"] == strain
        ].index

        tile_start = ds["peptide_metadata"]["nt_start"][strain_indices]
        tile_end = ds["peptide_metadata"]["nt_end"][strain_indices]
        enrichment = ds["counts"][sample][strain_indices]
        lines = [
            [(x1, y), (x2, y)] for x1, x2, y in zip(tile_start, tile_end, enrichment)
        ]
        lc = mc.LineCollection(lines, linewidths=1, label=strain, color=color_map[i])
        ax.add_collection(lc)

    ax.autoscale()
    sample_info = ds["sample_metadata"]["sample_info"][sample]
    seroepi_paper_id = ds["sample_metadata"]["seroepi_paper_id"][sample]
    notes = ds["sample_metadata"]["Notes"][sample]
    ax.set_title(
        f"sample_id : {sample},\n \
        sample_info : {sample_info},\n \
        seroepi_paper_id : {seroepi_paper_id},\n \
        Notes : {notes}"
    )
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_ylabel(f"Standardized fold enrichment")
    ax.set_xlabel(f"peptide tile, in order (genome position)")
    # ax.grid()
    plt.tight_layout()
    fig.savefig(f"{out}")
    plt.close(fig)
