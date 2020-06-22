# Phippery

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery)

## What is this?

**In short,**

[Phage Immunoprecipitation Sequencing](https://www.nature.com/articles/s41596-018-0025-6)
(PhIP-Seq)
Is a powerful protocol for investigating antibody binding specificities and potential pathogen epitopes.
This package provides a few tools to help collect data coming from a Nextflow pipeline seen
[here](https://github.com/matsengrp/phip-flow).
as well as source code to query the resulting
[xarray](http://xarray.pydata.org/en/stable/)
dataset.


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
Currently, the CLI is mostly used for as a utility to

```
Usage: phippery [OPTIONS] COMMAND [ARGS]...

  Some tools for PhIP-Seq data analysis. For help with a specific command,
  type phippery COMMAND -h

Options:
  -h, --help  Show this message and exit.

Commands:
  collect-phip-data    Collect sample and peptide metadata tables along
                       with...

  peptide-md-to-fasta  convert peptide metadata to fasta format.
```

