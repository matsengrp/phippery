# Phippery

<p>
  <img src="data/cartoons/Xarray_function.png" width="250">
</p>

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docker Repository on Quay](https://quay.io/repository/hdc-workflows/phippery/status "Docker Repository on Quay")](https://quay.io/repository/hdc-workflows/phippery)
[![build and test](https://github.com/matsengrp/phippery/actions/workflows/build-and-test.yaml/badge.svg?branch=main)](https://github.com/matsengrp/phippery/actions/workflows/build-and-test.yaml)
[![docs](https://github.com/matsengrp/phippery/actions/workflows/docs_pages_workflow.yml/badge.svg?branch=main)](https://github.com/matsengrp/phippery/actions/workflows/docs_pages_workflow.yml)
[![package](https://github.com/matsengrp/phippery/actions/workflows/package.yaml/badge.svg?branch=main)](https://github.com/matsengrp/phippery/actions/workflows/package.yaml)

A set of functions designed to query an 
[xarray DataSet](http://xarray.pydata.org/en/stable/user-guide/data-structures.html#dataset) 
object formatted to tie enrichment data with 
their respective row & column (peptide & sample) annotations. 

Please see the 
[documentation](https://matsengrp.github.io/phippery/) 
for further details.

## Developer Install

```
# using venv
python -m venv phippery_dev_env
source phippery_dev_env/bin/activate

# install
git clone https://github.com/matsengrp/phippery.git
(cd phippery && pip install -e ".[dev]")
```

