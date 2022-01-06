
.. note:: The software presented here is still under construction and 
    considered to be in the "Beta" stage of production. 
    Please expect and excuse innevitable changes, 
    for questions and/or suggestions, please feel welcome 
    to contact jgallowa (at) fredhutch.org


.. _sec_quick_start:

========
Examples
========

There are a few primary steps to PhIP-Seq analysis after the sequencing and
demultiplexing of samples. To address each of these, we provide

1.  A flexible `Nextflow automated pipeline <https://www.nextflow.io/>`_ 
    used for producing the 
    `raw enrichment data <TODO>`_ when provided 
    Next Generations Sequencing (demultiplexed `fastq files <TODO>`_) data, 
    as well as coupled `sample and peptide library annotation files <TODO>`_ 
    files, as input.

2.  A `Python <http://www.python.org/>`_ API/CLI for some useful query utilities.

3.  A `Streamlit <https://streamlit.io/>`_ interactive app for visualization 
    of aggregated enrichment of specified 
    sample-peptide binding events observed in a study.

Below, we'll start by running the pipeline on some example data.

.. _sec_soup_nutz:

Pan-CoV library example data
++++++++++++++++++++++++++++

.. below, we start with walk through a `example dataset <TODO>`_ which highlights the major
  steps involved in the workflow we typically suggest. 
  For brevity, the example does not go into much detail -- for more about each of the
  three steps described below, you may
  refer to the respective
  :ref:`Alignments Pipeline <sec_pipeline_intro>`,
  :ref:`Command Line Interface <sec_cli_intro>`, or
  :ref:`Interactive Visualization <sec_viz_intro>` pages.

The example dataset provided with this pipeline
is derived from a pre-COVID-19 heathy adult serum
sample, along with serum from a SARS-CoV-2 infected convalescent individual
both run in duplicate across two separate batches of Pan-Human CoV, full
proteome libraries. To read more on the dataset and respective input
inputs, see the more in-depth explanation in the 
`generalized alignment pipeline section <TODO>`_.
Huge thanks to
`Stoddard et al. 2021 <https://www.cell.com/cell-reports/fulltext/S2211-1247(21)00506-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124721005064%3Fshowall%3Dtrue>`_ and the 
`Overbaugh folks <TODO>`_.

.. _sec_align_soup_nutz:

Generating alignment counts data
++++++++++++++++++++++++++++++++

We provide some example pipeline input including the next generation
sequencing files for each of the six samples described
in the `sample table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/sample_table.csv>`_. We'll begin by aligning 
the fastq files to oligo library described by the 
Oligonucleotide encoding sequences in the 
`peptide table <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/peptide_table.csv>`_.

.. _sec_clone_template:

To run the pipeline example data, simply

.. code-block::

    nextflow run matsengrp/phip-flow --output_tall_csv true --output_wide_csv true -profile docker -resume

.. note::
    if you happen to get the message: *Project `matsengrp/phip-flow` currently is sticked on revision: <X> -- you need to specify explicitly a revision with the option `-r` to use it*. 
    Then run **nextyflow drop matsengrp/phip-flow**.
    See `this link <https://nf-co.re/usage/troubleshooting#warning-about-sticked-on-revision>`_ for more information about this error

Once the pipeline finishes, you'll see a summary of runtime stats.
Here we specified four parameters, two native to ``Nextflow`` 
(denoted with a single `-` prefix), and two which are specific to 
``phip-flow`` (double minus `--` symbols, for these).
the options ``output_tall_csv`` and ``output_wide_csv`` this options specifies one
of three optional output formats; tall csv, wide csv, and a pickle'd
binary xarray Dataset object. For more on these formats see this 
`great blog post <https://medium.com/w2hds/wide-tall-data-formats-423331ab5991>`_ 
on the topic.

This command ran the pipeline on the example dataset 
described above, the files can be viewed in the
`phip-flow git repo <https://github.com/matsengrp/phip-flow/tree/41_bin/data/pan-cov-example>`_.
In short, the workflow (1) used `bowtie` to align all the reads described in the 
sample annotation table, to the reference phage library described in the 
peptide table, (2) computed some usefule stats, and (3) formatted the data
into a single coherent dataset.
For more detail about the exact steps that were run, 
see the :ref:`nextflow pipeline page <>`.

Running on HPC (cluster)
++++++++++++++++++++++++

Above, we specified `-profile docker` which will assume you are running
this locally with *Docker* and *Nextflow* installed. 
for high performance computing systems, we can also specify
the `-profile cluster` option for running the default configurations
on a `slurm <https://slurm.schedmd.com/documentation.html>`_ cluster.
This option assumes the cluster has loaded modules or installs for 
*Singularity* and *Nextflow*. Here's an example script we might execute to run
the pipeline on the Fred Hutch Rhino machines:

.. code-block:: bash

    #!/bin/bash

    set -e
    source /app/lmod/lmod/init/profile

    module load nextflow
    module load Singularity
    export PATH=$SINGULARITYROOT/bin/:$PATH

    nextflow run matsengrp/phip-flow \
            --dump_tall_csv true \
            --dump_wide_csv true \
            --results "$(date -I)" \
            -profile cluster \
            -resume


Example results (tall csv)
++++++++++++++++++++++++++


Now, let's lets take a quick 
look at the results from the Pan-CoV example dataset that was run.
By default, the pipeline runs the Pan-CoV example data,
and writes the results out to a directory *results*.

::

  results
  ├── pickle_data
  │   └── data.phip
  ├── tall_data
  │   └── data-tall.csv
  └── wide_data
      ├── data_counts.csv
      ├── data_cpm.csv
      ├── data_enrichment.csv
      ├── data_peptide_annotation_table.csv
      ├── data_sample_annotation_table.csv
      └── data_size_factors.csv
  
  3 directories, 8 files
  
Let's take a look at how you might use **ggplot**
to visualize the data found in the tall formatted csv.
We'll start by plotting the individual sample enrichments, colored by
infection status.

.. code-block:: R

  library(ggplot2)
  library(dplyr)
  library(viridis)
  
  phip_data <- read.table(
          "results/tall_data/data-tall.csv", 
          header=TRUE, sep= ","
      ) %>%
      filter(Protein == "spike") %>%
      filter(Virus == "SARSCoV2") 
  
  # Plot
  p <- phip_data %>%
    ggplot(aes(
          x=Prot_Start, y=counts, 
          group=factor(sample_id), 
          color=factor(patient_status))
      ) +
      theme_bw() +
      geom_line() +
      ggtitle("Sars-CoV-2 Spike Protein Enrichments") +
      labs(y="# peptide alignments", color="infection status")
  ggsave(file="test.svg", plot=p, width=7, height=5)


.. figure:: images/example_counts.svg
  :width: 700
  :alt: example results
  :align: center

  Example data counts plotted as a function of location on the Spike
  protein of the 

Example results (wide csv)
++++++++++++++++++++++++++

Looking at the files in the wide format sub directory, we are given back the
peptide and sample annotation table's, both 
with an index (i.e. first) column "peptide_id" and "sample_id".
These indexes can simply be mapped back to the rows and columns
of each of the output enrichment matrices.
By default, the phip-flow pipeline outputs the raw counts, as well as
counts per million, and size factors (anders and huber, 2014 <TODO cite>)
normalizations of the matrix.
Let's use matplotlib's implot to plat the same samples as a heatmap.

.. code-block:: python3
  
    


.. _sec_cli_soup_nutz:

CLI for dataset query
+++++++++++++++++++++

With the binary dataset output (default)
and an installation of
the :ref:`phippery <sec_installation_phippery>`_ CLI tools,
we can run the some useful queries on the dataset to learn a little
about the dataset.

.. For all types of analysis outside of read alignment and visualization, 
  we recommend using the Command Line Interface (CLI) accessed using the 
  :program:`phippery` command.
  First, we'll take a look at the dataset using the 
  :program:`about` subcommand.

.. code-block::

  $ phippery about output/Pan-CoV-example.phip

The `about` will print information about 
the three primary aspects of a single dataset; Samples, Peptides, and Enrichment
Layers. For more about how the data is structured, 
see the :ref:`under the hood <sec_python_intro>`_ page.
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

Run the Visalization app
++++++++++++++++++++++++

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
