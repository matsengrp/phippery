.. _sec_introduction:

============
Introduction
============

.. note:: This documentation is incomplete and under development. If
    you would like to help, please open an issue or pull request at
    `GitHub <https://github.com/matsengrp/phippery/>`_.

This is the documentation for ``phippery``, a set of functions designed to help query
`PhIP-Seq <https://www.nature.com/articles/s41596-018-0025-6>`_ 
data in the form of an 
`xarray <http://xarray.pydata.org/en/stable/index.html>`_ 
DataSet which ties sample and 
peptide annotations to a enrichment matrix.



In a nutshell
-------------

We designed ``phippery`` to address a lack of standardized tools for analyzing counts the
*enrichment* tables obtained as a result of the protocol's sequence alignment heuristic.
There are a few primary steps to PhIP-Seq Data analysis post sequencing of the sample 
IP's: sample alignment to phage library, queries on the resulting enrichment matrix, and
data visualization to explore the data. 
Here, we have provided a 
`Nextflow <https://www.nextflow.io/>`_ pipeine, a
`Python <http://www.python.org/>`_ API/CLI, and a
`Streamlit <https://streamlit.io/>`_ visualization app 
which can be used separately, or in conjuction for the
rapid exploration of PhIP-Seq data.

TODO - Put full schematic, here.

First steps
-----------

Head to the :ref:`Installation <sec_installation>` page to get installation 
instructions for each of the three tools described above.

Motivation
----------

The primary data strucure resulting from PhIP-Seq experiments is an *enrichment matrix*, 
X, with i rows and j columns. 
Commonly, row index represents a peptide that is displayed on a phage,
and each column represents a sample that was mixed with the entire phage library. 
After sequencing and demultiplexing each sample, we align the reads to the 
oligonucleotide reference library to observe a
count of aligned reads to each peptide.

Outside of the enrichment matrix, each *sample* in an experiment as well as each *peptide*
in the phage library used have number of important annotations required when
performing analysis tasks like model fitting, normalizing, and differential selection.
Additionally, the comparison across groups of virus proteins and 
sample types is crucial in many experiments. For large sample size experiments, 
it can be difficult to cross reference each of these groups before and
after analysis. 
Here, we take advantage of the powerful 
`xarray <http://xarray.pydata.org/en/stable/>`_
approach to organizing all the Phip-Seq data along four primary coordinate 
dimensions which tie all sample/peptide enrichments to the respective annotations. 
Doing this allows us to store all the information without the error prone 
step of cross-checking separate dataframes, and without the
large storage scaling of using "Tall" dataframes.

Under the hood,

TODO - Explain more in detail how to use all three tools, probably a nice visual

