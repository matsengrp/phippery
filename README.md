# Phippery

A set of tools to organize, slice, and query data from PhIP-Seq style enrichment matrices in
a an effecient manner centered around the [xarray](http://xarray.pydata.org/en/stable/) infrastructure.

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery)
![build and test](https://github.com/matsengrp/phippery/workflows/build%20and%20test/badge.svg)

<p align="center">
  <img src="data/cartoons/Xarray_function.png" width="350">
</p>

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
cd phippery
conda env create -f environment.yml #follow prompts
conda activate phippery
```

install by
```
pip install .
```

Then run the tests:
```
pytest
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

