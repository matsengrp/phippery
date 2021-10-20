

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

In this page, we provide a few examples to highlight the features of each.

.. _sec_soup_nutz:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Pan-CoV library "soup-to-nuts"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below, we start with walk through a `example dataset <TODO>`_ which highlights the major
steps involved in the workflow we typically suggest. 
For brevity, the example does not go into much detail -- for more about each of the
three steps described below, you may
refer to the respective
:ref:`Alignments Pipeline <sec_pipeline_intro>`,
:ref:`Command Line Interface <sec_cli_intro>`, or
:ref:`Interactive Visualization <sec_viz_intro>` pages.

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


.. _sec_align_soup_nutz:

Step 1. Align fastq reads to peptide library
++++++++++++++++++++++++++++++++++++++++++++

We provide some example pipeline input including the next generation
sequencing files for each of the six samples described
in the `sample table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/sample_table.csv>`_. We'll begin by aligning 
the fastq files to oligo library described by the 
Oligonucleotide encoding sequences in the 
`peptide table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/peptide_table.csv>`_.

.. _sec_clone_template:

First, clone the `phip-flow-template <TODO>`_ 
(this may take a minute) and change the working directory to
`phip-flow-template/Pan-CoV-example-ds/`.

.. code-block:: bash

  » git clone git@github.com:matsengrp/phip-flow-template.git

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


Now that we have some input, 
we can use `Nextflow's git aware <TODO>`_ infrastructure to
run the bleeding edge script directly from the source 
`git repository <https://github.com/matsengrp/phip-flow>`_
on the Pan-CoV example files. The shell script inside the 
working directory, `run_phip_flow.sh` gives an example of
the command one might use to run the pipeline.

Running

.. code-block:: bash
  
  » /usr/bin/time ./run_phip_flow.sh

Produces the following output:

::

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

After the pipeline has completed it's run,
If using the default config file,
the output is an `xarray Dataset <TODO>`_
pickle dumped to a binary file `output/Pan-CoV-example.phip`

.. tip:: At anytime, this binary pickle dump
  can be converted to `tall`, or `wide` csv formats
  using the :program:`to-tall-csv` or 
  :program:`to-wide-csv` commands in the 
  :ref:`phippery CLI <sec_cli_intro>`
  For more information on the output data structure,
  see :ref:`under the hood <sec_python_intro>`.

.. _sec_cli_soup_nutz:

Step 2. CLI for dataset query
+++++++++++++++++++++++++++++

Once the alignments have run and we have our binary dataset files.
We can install and run the some queries on the dataset to learn a little
about the dataset.

For all types of analysis outside of read alignment and visualization, 
we recommend using the Command Line Interface (CLI) accessed using the 
:program:`phippery` command.
First, we'll take a look at the dataset using the 
:program:`about` subcommand.

.. code-block::

  $ phippery about output/Pan-CoV-example.phip

This will print information about the three primary aspects of a single dataset.
Primarily, it tells you what information is available in terms of the 
`Samples Table`,
`Peptide Table`,
and `Enrichment Layers`.

::

  Sample Table:
  -------------
  <class 'pandas.core.frame.DataFrame'>
  Int64Index: 6 entries, 124 to 540
  Data columns (total 10 columns):
   #   Column               Non-Null Count  Dtype
  ---  ------               --------------  -----
   0   seq_dir              6 non-null      string
   1   library_batch        6 non-null      string
   2   control_status       6 non-null      string
   3   participant_ID       4 non-null      string
   4   patient_status       4 non-null      string
   5   fastq_filename       6 non-null      string
   6   raw-total-sequences  6 non-null      Int64
   7   reads-mapped         6 non-null      Int64
   8   error-rate           6 non-null      Float64
   9   average-quality      6 non-null      Float64
  dtypes: Float64(2), Int64(2), string(6)
  memory usage: 552.0 bytes

Above we see our example dataset `sample table`. 
The information about
annotation feature data types, and missing information (NA) counts 
is provided by default.

As displayed, this dataset contains 6 samples, 
each with the annotations we fed to the pipeline
along with some alignment statistics.
While maybe not immediately useful, it's nice to know
which information you have available at any given time --
especially after we start slicing or grouping datasets. 

Further, you may want to know more detail about one of the annotation columns
at a time. The :program:`about-feature` will give you a useful description 
of the feature level distributions (categorical or numeric features), as well
as a few example queries for help indexing the dataset by this annotation feature
(we talk more about queries in the :ref:`next example <sec_neg_binom>`).
Let's take a look at our 
`reads mapped <http://www.htslib.org/doc/samtools-stats.html>`_ 
annotation feature:

::
  
  reads-mapped: Integer Feature:
  ---------------------------
  
  distribution of numerical feature:
  
  count         6.000000
  mean     359803.000000
  std      283811.764886
  min      122878.000000
  25%      147733.250000
  50%      234885.500000
  75%      597263.000000
  max      729431.000000
  Name: reads-mapped, dtype: float64
  
  
  Some example query statements:
  ------------------------------
  
  > "reads-mapped >= 359803"
  
  > "reads-mapped <= 359803"
  
  > "(reads-mapped >= 147733) and (reads-mapped <= 234885)"


.. Tip:: run `phippery -h` for a list of possible Commands. Additionally, you can run
    

.. _sec_viz_soup_nutz:

Step 3. Run the Visalization app
++++++++++++++++++++++++++++++++

Now that we have computed some normalizations on our samples, 
we can go ahead and use the binary dataset as input to the interactive
visualization app. Running the command below will open the all in your browser
search for all valid binary xarray formatted datasets in the current working
directory for visualization. 
Currently, you need to have the `.phip` binary file from above
in the same working directory where you cloned the streamlit app.
Navigate to the directory where you cloned
`phip-flow <https://github.com/matsengrp/phip-viz>`_ repository,
and createp a 
`symlink <https://kb.iu.edu/d/abbe>`_ 
to the binary `.phip` file. 

.. code-block::
    
    ln -s ../phip-flow-template/Pan-CoV-example-ds/output/Pan-CoV-example.phip ./
    streamlit run streamlit-app.py

The app will fire up your default (or most recently opened) browser
and you're ready to make your first visualizations!


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


