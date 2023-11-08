
.. _sec_pipeline_intro:

===================
Alignments Pipeline
===================

A flexible `Nextflow automated pipeline <https://www.nextflow.io/>`_ 
used for producing the 
:ref:`enrichment data <sec_pipeline_outputs>`
from PhIP-Seq data when provided the demultiplexed 
:ref:`fastq files <sec_input_fasta>`,
as well as annotation files for both the experimental
:ref:`samples <sec_sam_anno>` and 
:ref:`peptides <sec_pep_anno>` in the phage library being used.

.. figure:: images/dagt.svg
   :class: with-border
   :width: 600
   :alt: directed acyclic graph 
   :align: left

   **Directed Acyclic Graph (DAG)** 
   defining the execution order and dependencies of each individual
   processes. The major steps of the pipeline can be summarized as:
   (1) Build a ``Bowtie`` index from the relevant peptide oligos
   (2) Align each of the samples to the library reference using
   `Bowtie` end-to-end alignment allowing for up to N mismatches (default 2).
   The user specifies the size of both the reads and peptide,
   the low-quality end of the read are then trimmed to match
   the reference length before alignment.
   (3) Peptide counts for each sample alignment are obtained
   using ``samtools-idxstats`` (`Li et al., 2009 <https://doi.org/10.1093/bioinformatics/btp352>`_) in parallel
   to computing the common alignment stats such as
   raw total sequences, reads_mapped, error_rate, and average_quality, by default.
   (4) The resulting dataset containing the enrichment matrix,
   sample metadata, and peptide metadata are organized
   using the `xarray <https://xarray.pydata.org/en/stable/#>`_
   package (`Hamman and Hoyer, 2017 <http://doi.org/10.5334/jors.148>`_).
   (5) Optionally, summary statistics such as counts per million,
   size factors, fold enrichment, as well as model fits for estimates
   of significance are computed.
   (6) By default, the pipeline outputs all results
   computed above as a pickle dump binary of the xarray object
   with all combined results, as well as individual csv's for each statistic computed. 
   You may also export the results in ``PhIPData`` (Chen et al 2022) format
   via RDS file, and even tall CSV formats.

.. _sec_pipeline_inputs:

===========
Input files
===========

.. _sec_sam_anno:

Sample table 
++++++++++++

A CSV where one of the columns must be ``fastq_filename`` listing
all samples to be run through the pipeline.
By default, the pipeline assumes the reads are relative to
the project directory where the pipeline is being executed.

.. note:: If there is some other prefix for the filepaths,
    you may check out the ``--reads_prefix`` parameter.

As an example, let's assume there's some directory ``ngs/*`` containing all the
fastq files for a project. To organize these files (excluding barcode files) 
into a minimal sample table describing each of their relative paths, we might 
use the following command.

.. code-block:: bash
  
    (echo "fastq_filepath" && ls ngs/*R1*.gz)  > sample_table.csv

.. seealso:: An example of a simple sample table can be found 
    `here <https://github.com/matsengrp/phip-flow/blob/main/data/pan-cov-example/sample_table.csv>`_.

In addition to the file paths provided for each sample, 
you may include as many colorful annotations as you would
like so long as the CSV stays tidy. 
One caveat is that you should *not* include a column named ``sample_id``,
this feature name is reserved for integer id's
which are added to both the sample and peptide annotation tables
for internal and user organization of the annotations with 
enrichment tables.
Many of the ``phippery`` API utilities,
are generalized to help you index, and otherwise
manipulate the data to your liking using any combination
of these annotations, so go wild with annotations!

Internally, data types are handled through conversion to pandas data types
- so it's best to keep data types consistent
between the columns provided. For :ref:`missing data <sec_missing_data>`, 
we recommend empty strings, "", 
but "NaN" and "N/A" also work (hopefully) as expected.

.. note:: Some of the :ref:`optional workflows <sec_optional_workflows>`
    have additional required annotations, so keep an eye for those.

.. _sec_input_fasta:

.. note:: The fastq files pointed to by the sample table described above
    are assumed to have uniform (trimmed) read lengths.
    During alignment, this is enforced by reads being 
    trimmed on the 3' end to match the length specified 
    by the ``--oligo_tile_length`` parameter. 

See :ref:`pipeline parameters <sec_pipeline_params>` for more.

.. _sec_pep_anno:

Peptide table
+++++++++++++

A CSV where one of the columns must be "oligo" which
contains the oligonucleotide sequence encoding a peptide in
the phage library. 
Adapters are assumed to be encoded by lower case nucleotides
and are ultimately tossed from the output of the pipeline.
Conversely, the oligonucleotide encoding of the expressed peptide
should be upper case.
Similar to the sample annotation table, you may include any
annotations you like to the peptides (e.g. "Virus", "Strain", "Loci" etc)
*except* an annotation named ``peptide_id`` which is again reserved for
the pipeline execution.

.. seealso:: An example of a simple peptide table can be found 
    `here <https://github.com/matsengrp/phip-flow/blob/main/data/pan-cov-example/peptide_table.csv>`__.

.. _sec_pipeline_outputs:
  
================
Pipeline results
================

The primary use of this pipeline is to process raw sequencing data,
produce the peptide counts table, apply statistical methods 
(such as the :ref:`EdgeR <sec_edger>`), then combine and organize
the results from these workflows for the user to analyze however they wish.
By default the pipeline will produce the following outputs 

::

  results
  ├── pickle_data
  │   └── data.phip
  ├── rds_data
  │   └── PhIPData.rds
  └── wide_data
      ├── data_counts.csv.gz
      ├── data_cpm.csv.gz
      ├── data_edgeR_hits.csv.gz
      ├── data_edgeR_logfc.csv.gz
      ├── data_edgeR_logpval.csv.gz
      ├── data_peptide_annotation_table.csv.gz
      ├── data_sample_annotation_table.csv.gz
      └── data_size_factors.csv.gz

  4 directories, 11 files

see the :ref:`example page <sec_quick_start>` 
for a more detailed explanation of these outputs.

.. _sec_pipeline_params:

==========
Parameters
==========

Below, we describe each of the possible parameters that may be passed to the pipeline.
Parameters with a "*" next to the name must be provided values
explicitly in the ``nextflow run``, command unless 
you wish to be using the default values described below.
Otherwise, the parameter value is only required for relevant the 
:ref:`optional workflow <sec_optional_workflows>`.


``--sample_table``

- help: Table describing each input sample, minimally containing the column 'fastq_filepath' with the name of each file to be analyzed. Control samples are indicated with a value of 'beads_only' in the column 'control_status'.
- wb_type: file
- required: True

``--reads_prefix``

- help: Folder which contains the files listed in the sample table
- wb_type: folder
- required: True

``--read_length``

- help: Read length for alignment
- wb_type: integer
- default: 125

``--fastq_stream_func``

- help: Set this as 'cat' if fastq files not g'zipped
- wb_type: string
- default: zcat

``--peptide_table``

- help: Table describing each peptide in the library, minimally containing the column 'oligo' with the sequence used for each peptide
- wb_type: file
- required: True

``--peptide_tile_length``

- help: Peptide length for alignment
- wb_type: integer
- default: 117

``--dataset_prefix``

- help: String which is prepended to all output files
- wb_type: string
- default: data

``--output_pickle_xarray``

- help: Generate output files in xarray pickle format
- wb_type: bool
- default: True

``--output_tall_csv``

- help: Generate output files in tall CSV format
- wb_type: bool
- default: True

``--output_wide_csv``

- help: Generate output files in wide CSV format
- wb_type: bool
- default: True

``--n_mismatches``

- help: Number of mismatches allowed
- wb_type: integer
- default: 2

``--bowtie_optional_args``

- help: Other bowtie options
- wb_type: string
- default: --tryhard --nomaqround --norc --best --sam --quiet

``--replicate_sequence_counts``

- help: Flag for replicating counts for replicate sequences
- wb_type: bool
- default: True

.. _sec_optional_workflows:

===================
Optional Parameters
===================

We provide a popular (at least for us)
selection of the features found in the
:mod:`phippery` python API as optional during pipeline
execution. To run any one of these
optional workflows, you'll set the relevant
boolean flag parameter to true. 
Additionally, you may need to provide
certain annotation features
and factor levels in the sample or peptide
table.

Our `example pan-CoV dataset <https://github.com/matsengrp/phip-flow/tree/main/data/pan-cov-example>`__
includes library enrichment samples that
are appropriately annotated in the 
sample table, meaning we could
run the cpm enrichment workflow like so:

.. code-block:: bash

    (base) ubuntu phippery/phip-flow ‹V1.04*› » nextflow run main.nf -profile docker --run_cpm_enr_workflow true
    N E X T F L O W  ~  version 21.04.3
    Launching `main.nf` [distracted_banach] - revision: 9ea43df075
    P H I P - F L O W!
    Matsen, Overbaugh, and Minot Labs
    Fred Hutchinson CRC, Seattle WA
    ================================
    sample_table    : /home/jared/MatsenGroup/Projects/phippery/phip-flow/data/pan-cov-example/sample_table.csv
    peptide_table   : /home/jared/MatsenGroup/Projects/phippery/phip-flow/data/pan-cov-example/peptide_table.csv
    results         : /home/jared/MatsenGroup/Projects/phippery/phip-flow/results/


    executor >  local (29)
    [2c/30601d] process > ALIGN:validate_sample_table (1)    [100%] 1 of 1 ✔
    [1d/073399] process > ALIGN:validate_peptide_table (1)   [100%] 1 of 1 ✔
    [37/6937e7] process > ALIGN:generate_fasta_reference (1) [100%] 1 of 1 ✔
    [61/a636a9] process > ALIGN:generate_index (1)           [100%] 1 of 1 ✔
    [1d/454757] process > ALIGN:short_read_alignment (1)     [100%] 6 of 6 ✔
    [0d/320a4d] process > ALIGN:sam_to_counts (6)            [100%] 6 of 6 ✔
    [f4/687d71] process > ALIGN:sam_to_stats (6)             [100%] 6 of 6 ✔
    [c6/036a15] process > ALIGN:collect_phip_data (1)        [100%] 1 of 1 ✔
    [c0/3eb016] process > ALIGN:replicate_counts (1)         [100%] 1 of 1 ✔
    [e8/fae63a] process > STATS:counts_per_million (1)       [100%] 1 of 1 ✔
    [a9/81d0b7] process > STATS:size_factors (1)             [100%] 1 of 1 ✔
    [11/031d9e] process > STATS:cpm_fold_enrichment (1)      [100%] 1 of 1 ✔
    [-        ] process > STATS:fit_predict_neg_binom        -
    [-        ] process > STATS:fit_predict_zscore           -
    [7e/df19de] process > STATS:merge_binary_datasets        [100%] 1 of 1 ✔
    [c0/5b1faf] process > DSOUT:dump_binary                  [100%] 1 of 1 ✔
    [-        ] process > DSOUT:dump_wide_csv                -
    [-        ] process > DSOUT:dump_tall_csv                -
    [-        ] process > AGG:split_samples                  -
    [-        ] process > AGG:aggregate_organisms            -
    [-        ] process > AGG:join_organisms                 -

We can then use the :mod:`phippery.utils` to read in the data to take a look at the results.

.. code-block:: python

    >>> import phippery
    >>> ds = phippery.load("results/pickle_data/data.phip")
    >>> ds
    <xarray.Dataset>
    Dimensions:           (sample_id: 6, peptide_id: 10047, sample_metadata: 9,
                           peptide_metadata: 7)
    Coordinates:
      * sample_id         (sample_id) int64 0 1 2 3 4 5
      * peptide_id        (peptide_id) int64 0 1 2 3 4 ... 10043 10044 10045 10046
      * sample_metadata   (sample_metadata) object 'library_batch' ... 'average_q...'
      * peptide_metadata  (peptide_metadata) object 'Virus' ... 'Prot_Start'
    Data variables:
        counts            (peptide_id, sample_id) int64 0 1 0 3 0 3 ... 0 1 0 0 2 1
        sample_table      (sample_id, sample_metadata) object 'MEGSUB' ... 37.3
        peptide_table     (peptide_id, peptide_metadata) object 'OC43' ... 4052
        size_factors      (peptide_id, sample_id) float64 0.0 1.0 0.0 ... 2.29 1.0
        cpm               (peptide_id, sample_id) float64 0.0 94.85 ... 560.2 245.6
        enrichment        (peptide_id, sample_id) float64 0.02065 1.979 ... 5.093
    >>> ds.counts.to_pandas()
    sample_id   0  1  2  3  4  5
    peptide_id
    0           0  1  0  3  0  3
    1           7  2  3  0  5  1
    2           2  1  0  0  0  0
    3           0  2  0  0  0  0
    4           0  0  0  0  0  0
    ...        .. .. .. .. .. ..
    10042       0  1  0  0  0  0
    10043       1  2  0  0  0  0
    10044       0  1  0  0  0  0
    10045       8  0  0  0  1  0
    10046       0  1  0  0  2  1
    
    [10047 rows x 6 columns]
    >>> ds.enrichment.to_pandas()
    sample_id          0         1         2          3          4          5
    peptide_id
    0           0.020650  1.979353  0.020651  34.805577   0.020651  15.238485
    1           1.552418  0.447580  4.091275   0.002347   3.289535   0.578875
    2           1.328661  0.671337  0.007004   0.007004   0.007004   0.007004
    3           0.010433  1.989571  0.010433   0.010433   0.010433   0.010433
    4           1.000000  1.000000  1.000000   1.000000   1.000000   1.000000
    ...              ...       ...       ...        ...        ...        ...
    10042       0.020650  1.979353  0.020651   0.020651   0.020651   0.020651
    10043       0.666665  1.333336  0.006992   0.006992   0.006992   0.006992
    10044       0.020650  1.979353  0.020651   0.020651   0.020651   0.020651
    10045       1.997354  0.002643  0.002643   0.002643   0.742907   0.002643
    10046       0.020650  1.979353  0.020651   0.020651  11.589568   5.093262
    
    [10047 rows x 6 columns] 


BEER
++++

.. warning::
    This workflow has not been fully tested and may be very slow.
    For errors which may arise from the BEER workflow, we recommend
    that you direct questions to the BEER developers.
    If you would like to run BEER outside of the pipeline, note that
    by default the pipeline runs EdgeR and outputs
    those results into the 
    `PhIPData <https://bioconductor.org/packages/release/bioc/html/PhIPData.html>`__
    object file which can be directly loaded and used with the BEER library.

``--run_BEER``

- help: Flag for running edgeR and BEER using the infrastructure in the
  BEER pipeline. See the
  `R Vignettes <http://www.bioconductor.org/packages/release/bioc/vignettes/beer/inst/doc/beer.html>`_
  and `BEER Paper <https://academic.oup.com/bioinformatics/article/38/19/4647/6663763>`_
  for more on this method.
  Enrichments, EdgeR hits, and Annotations are
  tied into a `PhIPData <https://www.bioconductor.org/packages/release/bioc/html/PhIPData.html>`_
  object and exported to an RDS binary object file.
  The object file is then saved in the ``params.results`` directory
  below the ``rds_data/`` sub-directory.
  Additionally, these results will be tied back into the
  xarray object used by the Python phippery API,
  as well as any CSV outputs (wide & tall).
- wb_type: bool
- default: False


CPM Enrichment
++++++++++++++

``--run_cpm_enr_workflow``

- help: Flag for running the enrichment workflow using counts
    per million as a pre-processing step to fold enrichment
    of empricial IP samples over library abundance controls.
    If ``True``, we require that the sample annotation table
    provides a column "control_status" for which a subset of samples
    is labeled as "library" indicating the sample is a control
    For a description of how this function works in more detail,
    see :meth:`phippery.normalize.enrichment`.
- wb_type: bool
- default: False


Z-Score
+++++++

``--run_zscore_fit_predict``

- help: Flag for running Z-score enrichment analysis.
    This model fits to mock ip (bead only controls)
    and thus requires the sample annotation column "control_status"
    where mock IP\'s are marked "beads_only".
    Note that this method uses the 
    :func:`phippery.normalize.counts_per_million` function
    to normalize the data before fitting and estimating significance
    using :func:`phippery.modeling.zscore`.
    For more on this method, see 
    :ref:`the background modeling documentation <sec_background_modeling>`
- wb_type: bool
- default: False


VirScan Organism Summary
++++++++++++++++++++++++

This workflow will summarize the hits to epitopes (peptides)
across the groups across the proteome specified in the peptide annotation
table input to the pipeline.

Note that this analysis workflow was created with the Virscan
peptide library in mind, but could be used for any peptide assay being
analyzed. For example, you could run this workflow on the example data like so:

.. code-block:: bash

    nextflow run matsengrp/phip-flow -r V1.12 \
            -profile docker \
            --summarize_by_organism true \
            --peptide_seq_col "Prot" \
            --peptide_org_col "Virus" \
            --results "$(date -I)"

``--summarize_by_organism``

- help: Flag used to control the summary of results by organism
    Requires that the peptide table includes information regarding
    the source organism for each epitope. It is possible to annotate
    an epitope as being contained in multiple organisms by including
    multiple lines with the same peptide.
    When this flag is enabled, an additional output table will be
    produced (``aggregated_data/organism.summary.csv.gz``) which summarizes
    the number of epitopes with Z-scores above the threshold
    (``--zscore_threshold``, described below) for each organism.
    The sequence of each peptide is taken into account to filter out
    overlapping peptide hits.
    For any pair of peptides which overlap by more than the allowed
    number of amino acids (``--max_overlap``), only the higher-scoring
    peptide (in terms of Z-score) will be retained.
    Overlaps between peptides are determined by exact k-mer matching.
    A peptide is marked as a 'hit' when it is above the threshold in
    all replicates of that sample. When it is only above the threshold
    in a subset of replicates, it is marked as 'discordant'.
    The Epitope Binding Score (EBS) is also calculated for each peptide
    as the mean Z-score across all replicates from the same sample.
    At the organism level, the max and mean EBS is reported.
    Finally, all of those results are reported for the subset of
    epitopes which are marked as 'public' (using ``--public_epitopes_csv``),
    which indicates that there is independent experimental evidence supporting
    the presence of binding antibodies in a human population.
- wb_type: bool
- default: False

``--peptide_org_col``

- help: Column in the peptide table indicating the organism for each peptide
- wb_type: string
- default: organism

``--peptide_seq_col``

- help: Column in the peptide table containing the peptide sequence (used to match against public epitopes provided with ``--public_epitopes_csv``)
- wb_type: string
- default: seq

``--max_overlap``

- help: Maximum allowed overlap between detected peptides
- wb_type: integer
- default: 7

``--zscore_threshold``

- help: Minimum Z-score threshold
- wb_type: float
- default: 2.5

``--sample_grouping_col``

- help: Column in the sample table used for mapping replicates to samples
- wb_type: string
- default:

``--public_epitopes_csv``

- help: Optional, a CSV containing public epitopes
- wb_type: file

``--public_epitopes_col``

- help: In the public epitopes CSV, the column containing the translated amino acid sequence
- wb_type: string
- default: peptide_translate

``--nxf_profile``

- help: Profile used for resource allocation (options: standard / docker / cluster)
- wb_env: PROFILE
- wb_type: string
- default: standard

