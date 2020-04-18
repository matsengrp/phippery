# Phippery

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## What is this?

[Phip-seq](https://www.nature.com/articles/s41596-018-0025-6)
is a powerful protocol for investigating potential anitbody targets
(epitopes). Phage immunoprecipitation sequencing works by synthesizing an
array of proteins -- each of which could potentially model the epitope of
pathogen -- by cloning genetic "tiles" into a phage vector of choice.
We then take samples of serum antibodies from patients who may have
experienced a full
adapative immune response and thus, may have produced antibodies able to
identify epitodes of pathogens of interest. Using magnetic beads, we extract
the synthetic phages which have had a protein interaction with the
paratopes found on the serum antibodies of recovered patients. This is the 'immunoprecipitation'.

After an immuno-precipitation, the phages caught on the magnetic beads are
amplified and barcoded. We then sequence the pahges to see which protein
sequences, or 'peptidome tiles' are enriched for in the
immunoprecipitation step. The pipeline consists demultiplexing samples,
then aligning the peptides to the reference library in order to
generate a count of each specific peptide for each sample. The final
result is in the form of a matrix of counts, X, where
we have a column (i) for each sample and a row for each peptide (j).

From this experiment, there are obvious sources from which we can expect noise
in the resulting raw data including;
IgG content of a sample,
"immuno-dominant" antibodies,
sequencing bias (GC-content?),
phage display,
etc.

To analyze these data matrices, we are looking for
a peptide which consitantly has high enrichment.
However, parsing the signal from the noise means we need to be
clever in the experiental setup (controls, replicates etc),
keep track of important metadata
for _both_ samples and amungst the noise,
and finally, fit complex models which can leverage controls
and replicates for parsing signal from noise. To do this, a few
methods have been proposed:

**normalized enrichment**
TODO

**poisson modeling**
TODO

**pep z-score**
TODO

`phippery` is a set of tools which will organize the data using the
powerful `xarray` package to tie together important metadata, replicates,
and samples all in one dataset to be easily queried.

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

## Reference



