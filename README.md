# Phippery

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## What is this?

[phage immunoprecipitation sequencing](https://www.nature.com/articles/s41596-018-0025-6)
or PhIP-Seq,
is a powerful protocol for investigating potential anitbody targets (epitopes).
In short, phip-seq looks for protein-protein interactions between
the human body's own antibodies with artificial pathogens (phage vectors) i.e.
the epitope displayed by a pathogen of interest.
The key assumption being that recovered patients who have
undegone a full adaptive responce to a disease we care about
will have high enough volumes of antibodies to accurately observe an
interaction with the pathogen causing the disease.

More concretely, Phip-seq works by
synthesizing an array of proteins
-- each of which could potentially model the epitope of pathogen --
by cloning genetic "tiles" (oligonucleotides encoding peptides)
into a phage vector of choice.
samples from a patient's serum antibodies are mixed with
magnetic beads and phages
to capture the enrichment of phages-antibody interaction.

After an immunoprecipitation, the phages-inserted oligos
caught on the magnetic beads are
prepped for sequencing, amplified and barcoded.
Using Next geneation sequencing,
researchers can then quantify the number of
each unique oligos representing a peptide
captured in the immunoprecipitation step using short read alignment.

The pipeline consists demultiplexing the samples,
then aligning the peptides to the reference library in order to
generate a count of each specific peptide and sample combination. The final
result is in the form of a matrix, X, where hueristically
there is a row, i, for each peptide
and a column, j, for each sample.

Then `X[j][i]` = # of reads for peptide j, from sample i

From this protocol, there are a variety of sources
from which we can expect add noise to out data
These are including but not limited to;
IgG content of a sample,
"immuno-dominant" antibodies,
sequencing bias (GC-content?),
amplification bias,
phage display.

To analyze these results, we are looking for
a peptide which consitantly has high enrichment.
However, parsing the signal from noise described above means we need to be
clever in the experiental setup (controls, replicates etc),
bookkeeping of important metadata
and finally, our modeling.
Analysis methods for this type of data is described in *Methods*

`phippery` is a set of tools which will organize the data using the
powerful `xarray` package to tie together important metadata, replicates,
and samples -- all in one dataset to be easily queried.
Given this common datastrucure, we can then compile set of tools which all
take the same input to produce (hopefully interesting) results.
This is what we hope to accomplish here.

## How do I install it?

To install the API and command-line scripts at the moment,
it suggested you clone the repository, create a conda
environment from `environment.yaml`, and run the tests to make
sure everything is working properly.

```
git clone https://github.com/matsengrp/phippery.git
conda env create -f environment.yaml #follow prompts
conda activate phippery
```

to run the tests:
```
pytest
```
at the top level directory

## CLI

## Input Data

## Examples

## Analysis Methods

**normalized enrichment**
TODO

**poisson modeling**
TODO

**pep z-score**
TODO

## References



