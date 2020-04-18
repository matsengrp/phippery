# Phippery

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## What is this?

[phage immunoprecipitation sequencing](https://www.nature.com/articles/s41596-018-0025-6)
or Phip-seq,
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

After an immunoprecipitation, the phages caught on the magnetic beads are
amplified and barcoded.
Using short read sequencing,
researchers can then quantify the number of
each unique phage inserterted oligos
captured in the immunoprecipitation step.

The pipeline consists demultiplexing the samples,
then aligning the peptides to the reference library in order to
generate a count of each specific peptide for each sample. The final
result is in the form of a matrix of counts, X, where hueristically
there is a row, i, for each peptide.
and a column, j, for each sample.

Then `X[j][i] = # of reads for peptide j, in sample i`

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



