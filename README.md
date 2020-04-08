# phippery

This repository is currently a sandbox for creating tools and pipeline infrastructure relating to 
[PhIP-seq](https://www.nature.com/articles/s41596-018-0025-6) 
data. 
To get a taste of the data we're dealing with and the stuff we're aiming to create, have a look at 
`notebooks/phip_analysis_python.ipynb`. 

## Installation 

If you would like to run some of the prelim analysis yourself,

1. `git clone https://github.com/matsengrp/phippery` && cd phippery
2. obtain `empirical_data.tar.gz`, extract it (tar -xvf empirical_data.tar.gz) in the top-level directory
3. create a new `conda` environment with the necessary dependencies `conda env create -f environment.yml`
4. then simply run jupyter `jupyter notebook` and open `notebooks/phip_analysis_python.ipynb`.




