#!/usr/bin/env python 
# import numpy as np

import pandas as pd
import argparse
import scipy.stats as st
import re
import sys
import os
from matplotlib import pyplot as plt

if __name__ == "__main__":

    # TODO Finish
    parser = argparse.ArgumentParser(
        description='Some python code \
        to run the PhiP-seq fold enrichment analysis. This requires \
        three different files, 1. the peptide counts file, peptide  \
        metadata, and sample metadata. Examples of each of these \
        file can be found ../test_files/.')

    parser.add_argument(
        '-counts', type=str,
        help='The csv counts matrix: The first column is the \
        peptide id to be matched in the peptide metadata first column. \
        The first row should indicate the sample number for each \
        technical replicate.')

    parser.add_argument(
        '-peptide_metadata', type=str,
        help='This should be a csv where the first column is the \
        id for each respective peptide sequence in the counts data.')

    parser.add_argument(
        '-sample_metadata', type=str,
        help='This csv file will contain the respective metadata \
        for each of our samples. The first column should match the \
        respective sample in the first row of the columns.')

    parser.add_argument(
        '-corr_th', type=float, default=0.8,
        help='To analyse the quality of each sample during lab \
        prep, we check the correlation between the samples technical \
        replicates. This argument is a float specifying how high \
        the correlation should be before throwing a sample out of \
        our analysis')

    parser.add_argument(
        '-out_dir', type=str, default="./fold_analysis_output",
        help='The directory you would like to put all analysis \
        figures and possible intermediate files.')

    args = parser.parse_args()
    if not os.path.exists(args.out_dir): os.mkdir(args.out_dir) 

    # check all input : aka match up all respective
    # metadata samples and peptide ids.
    counts = pd.read_csv(args.counts, header=0, index_col=0)
    peptide_metadata = pd.read_csv(args.peptide_metadata, header=0, index_col=0)

    # TODO, clean up the file and use this once we have all analysis.
    # sample_metadata = pd.read_csv(args.sample_metadata, header=0, index_col=0)
    
    # TODO just make a function which returns all counts and metadata 
    # after error checking the input
    # lets make sure we have the same indexing going on here.
    assert(set(counts.index) == set(peptide_metadata.index))

    # Now we need to look at correlation between the technical replicates
    # TODO there's certainly a cleaner way to do this.

    ## EXTRACT SAMPLE INFO () ##

    # TODO extract this from sample metadata.
    beads_only = "35"

    library_control = [re.match("\d+\.input",f"{sample_rep}")[0] 
        for sample_rep in counts.columns if re.match("\d+\.input",f"{sample_rep}") != None][0]

    print(library_control)

    technical_rep_1 = [re.match("\d+\.[1]",f"{sample_rep}")[0] 
        for sample_rep in counts.columns if re.match("\d+\.[1]",f"{sample_rep}") != None]
    technical_rep_2 = [re.match("\d+\.[2]",f"{sample_rep}")[0] 
        for sample_rep in counts.columns if re.match("\d+\.[2]",f"{sample_rep}") != None]

    print(f"technical replicates 1 are: {technical_rep_1}")
    print(f"technical replicates 2 are: {technical_rep_2}")

    ### TRANSFORM COUNTS - EXTRACT TECH REP CORR () ###

    correlations, samples = [], []
    for i, (t1, t2) in enumerate(zip(technical_rep_1, technical_rep_2)):
        
        # let's make sure we're looking at the same sample
        assert(t1.split(".")[0] == t2.split(".")[0])
        sample = t1.split(".")[0]

        # Now, lets determine which samples to drop,
        # and which samples we can sum replicates over
        # based off of correlation.
        corr = st.pearsonr(counts[t1], counts[t2])
        correlations.append(corr)
        samples.append(sample)
            
        # if the technical replicates are correlated enough then
        # sum the counts, add 10, and insert new column.
        print(f"sample {sample} has technical rep correlation: {corr[0]} \
            with a p-value of {round(corr[1],5)}")
        if (corr[0] > args.corr_th and corr[1] <= .05) or sample == beads_only:
            counts[sample] = counts[t1] + counts[t2] + 20
        else:
            print(f"dropping sample {sample}, b/c corr was: {corr[0]}")

        counts.drop([t1,t2], axis=1, inplace=True)

    fig, ax = plt.subplots(1)
    ax.bar(samples, [c[0] for c in correlations])
    ax.axhline(y=args.corr_th)
    fig.set_size_inches(10,5)
    fig.savefig(os.path.join(args.out_dir,"technical_replicates_correlation.pdf"))

    ### NORMALIZE ENRICHMENTS () ###
    """
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
    print(f"library control is: {library_control} where it should be 37.input" )

    print(f"column sums :\n {counts.sum(axis=0)}")
    #counts = counts.astype("float64")
    frequencies = counts.truediv(counts.sum(axis=0), axis=1)
    print(f"freq:\n {frequencies.head()}")
    enrichment = frequencies.truediv(frequencies[library_control], axis=0)
    print(f"enrichment:\n {enrichment.head()}")
    standardized_enrichment = enrichment.subtract(enrichment[beads_only], axis=0)
    print(f"standardized enrich:\n {standardized_enrichment.head()}")

    ### Extract ZIKA, DENGUE, and HIV specific strains for "sample".
    sample = "10"
    strain1_HIV = "HIV_Env_BF520_W14_C2"
    strain2_HIV = "HIV_Env_BG505_W6_C2"

    strain_1_meta = peptide_metadata[peptide_metadata["Virus_Strain"] == strain1_HIV]
    strain_2_meta = peptide_metadata[peptide_metadata["Virus_Strain"] == strain2_HIV]

    for sample in standardized_enrichment.columns:
        print(sample)
        #print(strain_1_meta.index)
        #print(strain_1_meta.shape[0])
        print(standardized_enrichment[sample][strain_1_meta.index])
        #sys.exit()
        fig, ax = plt.subplots()
        ax.plot(list(range(strain_1_meta.shape[0])),standardized_enrichment[sample][strain_1_meta.index])
        ax.plot(list(range(strain_2_meta.shape[0])),standardized_enrichment[sample][strain_2_meta.index])
        plt.show()
        #print(strain_2_meta)

        
        
        

        
        
        
        
        




















