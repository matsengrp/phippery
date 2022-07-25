.. _sec_cli_intro:

======================
Command Line Interface
======================

The CLI is written using the
`Click <https://click.palletsprojects.com/en/8.0.x/>`_
library, and thus both `phippery -h`, and `phippery COMMAND -h` will provide
the same information provided below

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

The **about** command will print information about 
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


.. click:: cli:cli
   :prog: phippery
   :nested: full
