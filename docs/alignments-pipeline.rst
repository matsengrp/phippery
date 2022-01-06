
.. _sec_pipeline_intro:

===================
Alignments Pipeline
===================

A flexible `Nextflow automated pipeline <https://www.nextflow.io/>`_ 
used for producing the 
`raw enrichment data <TODO>`_ when provided 
Next Generations Sequencing (demultiplexed `fastq files <TODO>`_) data, 
as well as coupled `sample and peptide library annotation files <TODO>`_ 
files, as input.


.. image:: images/dag.svg
  :width: 400
  :alt: Alternative text
  :align: center

Above, 

.. _sec_pipeline_inputs:

Here, we follow the same
:ref:`soup-to-nuts example (step 1) <sec_align_soup_nutz>` 
with a little more detail. We'll assume you've already
cloned the template repository.

.. code-block:: bash

  » git clone git@github.com:matsengrp/phip-flow-template.git

.. tip:: If you would like to retain a copy of the Nextflow 
  script locally for modification, use the `--recurse-submodules` flag.

++++++
Inputs
++++++

In the `Pan-CoV-example-ds/` directory we see a few files which define the complete input
into the alignment pipeline. 

.. code-block:: bash

  (base) ubuntu phippery/sandbox » cd phip-flow-template/Pan-CoV-example-ds
  (base) ubuntu phip-flow-template/Pan-CoV-example-ds ‹main› » tree -L 1
  .
  ├── NGS
  ├── peptide_table.csv
  ├── phipflow_docker.config
  ├── run_phip_flow.sh
  └── sample_table.csv

  1 directory, 4 files

The sample and peptide table are there to define crucial information
about individual, sequenced and demultiplexed sample IP's (or controls),
and individual peptides in the phage display library being used, respectively.

In the sample table, we have a unique integer id for each sample being aligned. 
In addition to the id, wemust have, at minimum,
annotations (data columns) which direct the `Nextflow` 
pipeline to each of the individual demultiplexed fastq files. 

In this example, we have pointed each of the sample file in `NGS/` directory.
Currently, we do this by combining `seq_dir` and `fastq_filename` columns.

.. code-block::

  » cat sample_table.csv | cut -d "," -f 1,2,7

::

  sample_id,seq_dir,fastq_filename
  540,NGS/20-05-15-cov2-ex4b/,4B-rep1-27-library_S26_L001_R1_001_sub.fastq.gz
  490,NGS/20-05-14-cov2-ex4a/,4A-rep1-27-library_S27_L001_R1_001_sub.fastq.gz
  124,NGS/20-05-27-cov2-ex5a/,rep1-15_S15_L001_R1_001.fastq.gz
  218,NGS/20-06-02-cov2-ex5b/,ex5b-rep1-15_S15_L001_R1_001.fastq.gz
  525,NGS/20-05-15-cov2-ex4b/,4B-rep1-18_S18_L001_R1_001.fastq.gz
  472,NGS/20-05-14-cov2-ex4a/,4A-rep2-18_S45_L001_R1_001.fastq.gz

.. warning::
  The sample_id's are always the first column in a sample table, and remain unique
  integers of your choosing when creating your dataset. :program:`phippery` 
  will maintain the integrity of these 
  identifiers throughout all analysis. 
  However, they will always be sorted when organized into the binary xarray 
  structure. The same is true for peptide id's


We then make sure that the filepaths above match the file structure 
of our NGS data. 

.. code-block::

  NGS
  ├── 20-05-14-cov2-ex4a
  │   ├── 4A-rep1-27-library_S27_L001_R1_001_sub.fastq.gz
  │   └── 4A-rep2-18_S45_L001_R1_001.fastq.gz
  ├── 20-05-15-cov2-ex4b
  │   ├── 4B-rep1-18_S18_L001_R1_001.fastq.gz
  │   └── 4B-rep1-27-library_S26_L001_R1_001_sub.fastq.gz
  ├── 20-05-27-cov2-ex5a
  │   └── rep1-15_S15_L001_R1_001.fastq.gz
  └── 20-06-02-cov2-ex5b
      └── ex5b-rep1-15_S15_L001_R1_001.fastq.gz

      4 directories, 6 files


.. tip:: For organzing fastq files that may be scattered among alarge file sysytem,
    Nextflow will follow `symbolic links <https://kb.iu.edu/d/abbe>`_ 
    pointed at by the Sample Table.

.. tip:: the file 
  `phip-flow-template/Pan-CoV-example-ds/phipflow_docker.config`
  contains all the relevent settings for running the alignment 
  pipeline using only the installs described above on any sufficient
  laptop. For more custom settings,
  see the `Nextlfow configuration documentation 
  <https://www.nextflow.io/docs/latest/config.html#configuration>`_.



.. code-block:: bash

  #!/bin/bash
  set -e
  
  /usr/bin/time nextflow  \
    -C phipflow_docker.config \
    run matsengrp/phip-flow/PhIP-Flow.nf -r main \
    -with-report ./output/nextflow_report.html \
    -work-dir ./output/work/ \
    -resume



--------------------------------
Next Generation Sequencing files
--------------------------------

------------------
Configuration File
------------------

---------------------
Input File Validation
---------------------

.. _sec_pipeline_outputs:

+++++++
Outputs
+++++++

..
  .. _sec_neg_binom:
  
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Fitting a Negative Binomial model to mock IP's
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  Coming soon ...
  
  
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  Calculating Fold enrichment with library 
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  Coming soon ...
  
  
  ^^^^^^^^^^^^^^^^^^^^^^
  Differential Selection
  ^^^^^^^^^^^^^^^^^^^^^^
  
  Coming soon ...


