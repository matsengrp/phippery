# Phippery

## THIS PACKAGE IS CURRENTLY UNDER ACTIVE DEVELOPMENT. 
**We expect a more stable version sometime in January, 2022. Stay tuned!**

<p>
  <img src="data/cartoons/Xarray_function.png" width="250">
</p>

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Docker Repository on Quay](https://quay.io/repository/matsengrp/phippery/status "Docker Repository on Quay")](https://quay.io/repository/matsengrp/phippery)
[![build and test](https://github.com/matsengrp/phippery/workflows/build%20and%20test/badge.svg)](https://github.com/matsengrp/phippery/blob/master/.github/workflows/build-and-test.yaml)

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

