

.. _sec_quick_start:

===========
Quick Start
===========

After installation of the requirements described above, you're all set to run
the entire analysis workflow from fastq files, and annotation tables, to running
the interactive exploration app. 
Below, we run the entire pipeline on an
example dataset on the COVID-19 humoral responce, subsetted from 
`Stoddard et al. 2021 <https://www.cell.com/cell-reports/fulltext/S2211-1247(21)00506-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124721005064%3Fshowall%3Dtrue>`_. 


.. note:: Huge thanks to
    `Overbaugh Lab <https://research.fredhutch.org/overbaugh/en.html>`_
    for data and feedback.

TODO - Tree

Step 1. Align Fastq reads to peptide library
++++++++++++++++++++++++++++++++++++++++++++

TODO - Finish

Please be sure the phip-flow repository has been cloned and your
current working directory is the top-level directory of that repo. 
Next, run the nextflow pipeline locally using docker containers and native cpu

::

  $ nextflow run matsengrp/phip-flow/PhIP-Flow.nf -C foo-phip.config -o bar.phip

The output HDF5 formatted 

Step 2. CLI for dataset query
+++++++++++++++++++++++++++++

TODO - Finish

In the environment where ``phippery`` and ``phipviz`` have been installed,
we may start by summarizing the dataset that was just produced by the pipeline above

::

  $ phippery summarize bar.phip


Next,

.. note:: For a list of other 
    commands available by the click CLI:  `$ phippery -h`  

Step 3. Run the Visalization app
++++++++++++++++++++++++++++++++

TODO - Finish

..
  https://www.section.io/engineering-education/how-to-deploy-streamlit-app-with-docker/

Run the phip-viz app, using the normalized sample data

::

  $ stream run ...

..
  Docker
  ''''''
  ::
  
    $ docker pull quay.io/matsengrp/phip-viz
    $ docker run -it -v bar.phip:/ quay.io/matsengrp/phip-viz:latest
