

.. _sec_quick_start:

========
Examples
========

.. note:: The software presented here is still under construction and 
    considered to be in the "Beta" stage of production. 
    Please expect and excuse innevitable changes, 
    for questions and/or suggestions, please feel welcome 
    to contact jgallowa (at) fredhutch.org

There are a few primary steps to PhIP-Seq data analysis after sequencing of the 
demultiplexing the samples. To address each of these, we provide

1.  A flexible `Nextflow automated pipeline <https://www.nextflow.io/>`_ 
    used for producing the 
    `raw enrichment data <TODO>`_ when provided 
    Next Generations Sequencing (demultiplexed `fastq files <TODO>`_) data, 
    as well as coupled `sample and peptide library annotation files <TODO>`_ 
    files, as input.

2.  A `Python <http://www.python.org/>`_ API/CLI for many useful queries, 
    normalizations, and model fitting procedures

3.  A `Streamlit <https://streamlit.io/>`_ interactive app for visualization 
    of aggregated enrichment of specified 
    sample-peptide binding events observed in a study.


Below, we start with walk through a `example dataset <TODO>`_ which highlights the major
steps involved in the workflow we typically suggest. 
For brevity, the example is quite minimal -- If you would like more detail
on any of the individual three steps described below, you may
refer to the respective
:ref:`Alignments Pipeline <sec_pipeline_intro>`,
:ref:`Command Line Interface <sec_cli_intro>`, or
:ref:`Interactive Visualization <sec_viz_intro>` section, as necessary.

.. note::
  The example dataset is derived from a pre-COVID-19 heathy adult serum
  sample, along with serum from a SARS-CoV-2 infected convalescent individual
  both run in duplicate across two separate batches of Pan-Human CoV, full
  proteome libraries. To read more on the dataset and respective input
  inputs, see the more in-depth explanation in the 
  `generalized alignment pipeline section <TODO>`_.
  Huge thanks to
  `Stoddard et al. 2021 <https://www.cell.com/cell-reports/fulltext/S2211-1247(21)00506-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124721005064%3Fshowall%3Dtrue>`_ and the 
  `Overbaugh folks <TODO>`_.

.. _sec_soup_nutz:

Pan-CoV library "soup-to-nuts"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _sec_align_soup_nutz:

Step 1. Align fastq reads to peptide library
++++++++++++++++++++++++++++++++++++++++++++

The example files we provide include next generation
sequencing files for each of the six samples described
in the `sample table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/sample_table.csv>`_. We'll begin by aligning 
the fastq files to oligo library described by the 
Oligonucleotide encoding sequences in the 
`peptide table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/peptide_table.csv>`_.

.. _sec_clone_template:

Next, clone the `phip-flow-template <TODO>`_  and change the working directory to
`phip-flow-template/Pan-CoV-example-ds/`

.. code-block:: bash

  (base) ubuntu phippery/sandbox » /usr/bin/time git clone git@github.com:matsengrp/phip-flow-template.git
  Cloning into 'phip-flow-template'...
  remote: Enumerating objects: 107, done.
  remote: Counting objects: 100% (107/107), done.
  remote: Compressing objects: 100% (76/76), done.
  remote: Total 107 (delta 51), reused 85 (delta 29), pack-reused 0
  Receiving objects: 100% (107/107), 182.88 MiB | 16.02 MiB/s, done.
  Resolving deltas: 100% (51/51), done.
  1.38user 3.10system 0:21.38elapsed 21%CPU (0avgtext+0avgdata 243572maxresident)k
  8inputs+753744outputs (0major+88045minor)pagefaults 0swaps

.. code-block:: bash

  (base) ubuntu phippery/sandbox » cd phip-flow-template/Pan-CoV-example-ds

.. tip:: the file 
  `phip-flow-template/Pan-CoV-example-ds/phipflow_docker.config`
  contains all the relevent settings for running the alignment 
  pipeline using only the installs described above on any sufficient
  laptop. For more custom settings,
  see the `Nextlfow configuration documentation 
  <https://www.nextflow.io/docs/latest/config.html#configuration>`_.


Finally, we can use `Nextflow's git aware <TODO>`_ infrastructure to
run the bleeding edge script directly from the source 
`git repository <https://github.com/matsengrp/phip-flow>`_
on the Pan-CoV example files. The shell script inside the 
working directory, `run_phip_flow.sh` gives an example of
the command one might use to run the pipeline.

.. code-block:: bash

  #!/bin/bash
  set -e
  
  /usr/bin/time nextflow  \
    -C phipflow_docker.config \
    run matsengrp/phip-flow/PhIP-Flow.nf \
    -with-report ./output/nextflow_report.html \
    -work-dir ./output/work/ \
    -resume
    
.. tip:: If you would like to retain a copy of the Nextflow 
  script locally for modification, use the `--recurse-submodules` flag.

.. code-block:: bash
  
  (phippery) ubuntu phip-flow-template/Pan-CoV-example-ds ‹main*› » /usr/bin/time ./run_phip_flow.sh
  N E X T F L O W  ~  version 20.04.1
  Pulling matsengrp/phip-flow ...
  Already-up-to-date
  Launching `matsengrp/phip-flow` [lethal_brown] - revision: 1dfb4d69a8 [master]
  executor >  local (21)
  [be/cfbe3a] process > generate_fasta_reference (1) [100%] 1 of 1 ✔
  [1c/e5bc1b] process > generate_index (1)           [100%] 1 of 1 ✔
  [a7/b40db3] process > short_read_alignment (2)     [100%] 6 of 6 ✔
  [47/8bd2e8] process > sam_to_stats (6)             [100%] 6 of 6 ✔
  [f4/5512e5] process > sam_to_counts (6)            [100%] 6 of 6 ✔
  [4d/e73232] process > collect_phip_data (1)        [100%] 1 of 1 ✔
  Completed at: 12-Oct-2021 02:35:14
  Duration    : 1m 2s
  CPU hours   : (a few seconds)
  Succeeded   : 21

The output is an `xarray Dataset <TODO>`_
formatted using `net CDF <TODO>`_.
For more information on the output data structure,
see :ref:`under the hood <sec_python_intro>`.

.. _sec_cli_soup_nutz:

Step 2. CLI for dataset query
+++++++++++++++++++++++++++++

Once the alignments have run and we have our binary dataset files.
We can install and run the some queries on the dataset to learn a little
about the dataset.


.. _sec_viz_soup_nutz:

Step 3. Run the Visalization app
++++++++++++++++++++++++++++++++

Now that we ave computed some normalizations of interest on our samples, 
we can go ahead and use the binary net CDF dataset as input to the interactive
visualization app
