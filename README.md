# Phippery

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## What is this?

**In short,**

[Phage Immunoprecipitation Sequencing](https://www.nature.com/articles/s41596-018-0025-6)
(PhIP-Seq)
is a powerful protocol for investigating potential antibody targets (epitopes).
PhIP-Seq looks for protein-protein interactions between
the human body's own antibodies with artificial pathogens (phage vectors) i.e.
the epitope displayed by a pathogen of interest.
The key assumption being that recovered patients who have
undergone a full adaptive response to a disease we care about
will have high enough volumes of antibodies to accurately observe an
interaction with the pathogen causing the disease. Here, we provide tools
to store the resulting data,
as well as compile a set of tools which take
this as input and produce results.

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

To run the tests:
```
pytest
```
at the top level directory.

If tests pass, you can install by:
```
python setup.py install
```
or, in developer mode:
```
pip install -e .
```

## CLI

`phippery` uses
[click](https://click.palletsprojects.com/en/7.x/) as a CLI manager. This means
that phippery has nested help pages for each command available.

```
(phippery) jaredgalloway@FortDellis:~/projects/phip/fork/phippery$ phippery -h
Usage: phippery [OPTIONS] COMMAND [ARGS]...

  Some tools for PhIP-Seq data analysis. For help with a specific command,
  type phippery COMMAND -h

Options:
  -h, --help  Show this message and exit.

Commands:
  collect-phip-data  Collect and merge counts/metadata into one dataset to be
                     pickle dump'd

  fold-analysis      perform standardized fold enrichment analysis on phip
                     dataset
```

## API

FINISH

Under the hood, `phippery` cli is calling a set of functions designed to
either collect raw counts data and metadata to structure into a `phip_dataset`
(`collect-phip-data` command), or run analysis onn that same dataset object.
This dataset is really just a dictionary (or xarray Dataset Object, if you desire)
with keys containing pandas DataFrames for "counts", "sample_metadata", and
"peptide_metadata".

## Input Data

FINISH

There are three main inputs `phippery` need to do create a `phip_dataset` to be used
by all analysis tools; a directory with each counts for each sample in the
experiment,  

**counts**

**sample metadata**

**peptide metadata**

## Examples

collect then compute fold analysis

```
set -eux
EMP_DIR=zika_denv_hiv_empirical
COUNTS=${EMP_DIR}/counts/
S_META=${EMP_DIR}/sample_metadata.tsv
P_META=${EMP_DIR}/peptide_metadata.tsv

PHIP_DS=phip_ds.phip
PHIP_DS_SE=phip_ds_se.phip

phippery collect-phip-data -c ${COUNTS} -s_meta ${S_META} -p_meta ${P_META} -o ${PHIP_DS}  
phippery fold-analysis -d ${PHIP_DS} -mock 35 -lib 37 -o ${PHIP_DS_SE}
rm ${PHIP_DS}

```

this results in a file `phip_ds_se.phip` which can be read and plotted using
the phippery API:

```python
import phippery.utils as ut
import pickle as pk
import pandas as pd

ds = pk.load(open("phip_ds_se.phip","rb"))

ut.plot_peptide_enrichment_by_nt_position(
    ds=ds,
    strain_pattern="HIV.+",
    sample=1,
    out=f"plots/test.pdf",
    cmap="magma"
)
```

which will result in a plot like this:

![alt text](<data/plots/example.png>)

## Analysis Methods

**normalized enrichment**
TODO

**poisson modeling**
TODO

**pep z-score**
TODO

## More PhIP-seq detail

PhIP-Seq works by
synthesizing an array of proteins
-- each of which could potentially model the epitope of pathogen --
by cloning genetic "tiles" (oligonucleotides encoding peptides)
into a phage vector of choice.
Samples from a patient's serum antibodies are mixed with
magnetic beads and phages
to capture the enrichment of phages-antibody interaction.

After an immunoprecipitation, the phages-inserted oligos
caught on the magnetic beads are
prepped for sequencing, amplified and barcoded.
Using Next generation sequencing,
researchers can then quantify the number of
each unique oligos representing a peptide
captured in the immunoprecipitation step using short read alignment.

The pipeline consists demultiplexing the samples,
then aligning the peptides to the reference library in order to
generate a count of each specific peptide and sample combination. The final
result is in the form of a matrix, X, where heuristically
there is a row, i, for each peptide
and a column, j, for each sample.

Then `X[j][i]` = # of reads for peptide j, from sample i

From this protocol, there are a variety of sources
from which one can expect add noise to the resulting data
These are including but not limited to;
IgG content of a sample,
"immuno-dominant" antibodies,
sequencing bias (GC-content?),
amplification bias, and
phage display (total amount of each phage in the phage library).

To analyze these results, we are looking for
a peptide which consistently has high enrichment.
However, parsing the signal from noise described above means we need to be
clever in the experimental setup (controls, replicates etc),
bookkeeping of important metadata
and finally, our modeling.
Analysis methods for this type of data is described in **Methods**

`phippery` is a set of tools which will organize the data using the
powerful `xarray` package to tie together important metadata, replicates,
and samples -- all in one dataset to be easily queried.
Given this common data structure, we can then compile set of tools which all
take the same input to produce (hopefully interesting) results.
This is what we hope to accomplish here.


## References

TODO
1. nature protocols
2. uri phip-stap
3. analysis as we go ...



