
.. _sec_pipeline_intro:

===================
Alignments Pipeline
===================

A flexible `Nextflow automated pipeline <https://www.nextflow.io/>`_ 
used for producing the 
:ref:`enrichment data <sec_pipeline_outputs>`
from phip-seq data when provided the demultiplexed 
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
   using ``samtools-idxstats`` (Li et al., 2009) in parallel
   to computing the common alignment stats such as
   raw total sequences, reads_mapped, error_rate, and average_quality, by default.
   (4) The resulting dataset containing the enrichment matrix,
   sample metadata, and peptide metadata are organized
   using the `xarray <https://xarray.pydata.org/en/stable/#>`_
   package (Hamman and Hoyer, 2017).
   (5) Optionally, summary statistics such as counts per million,
   size factors, fold enrichment, as well as model fits for estimates
   of significance are computed.
   (6) By default, the pipeline outputs all results
   computed above as a pickle dump binary of the xarray object
   described. Optionally, you may export the results in either wide, or
   tall CSV formats.

.. _sec_pipeline_inputs:

===========
Input files
===========

.. _sec_sam_anno:

Sample table 
++++++++++++

A CSV where one of the columns must be "fastq_filename" listing
all samples to be run through the pipeline.
By default, the pipeline assumes the reads are relative to
the project directory where the pipeline is being executed.

.. note:: If there is some other prefix for the filepaths,
    you may check out the ``--reads_prefix`` parameter.

As an example, let's assume there's some directory *ngs/* containing all the
fastq files for a project. To organize these files (excluding barcode files) 
into a minimal sample table describing each of their relative paths, we might 
use the following command.

.. code-block:: bash
  
    (echo "fastq_filepath" && ls ngs/*R1*.gz)  > sample_table.csv

.. seealso:: An example of a simple sample table can be found 
    `here <https://github.com/matsengrp/phip-flow-template/blob/main/Pan-CoV-example-ds/sample_table.csv>`_.

In addition to the filepaths provided for each sample, 
you may include as many colorful annotations as you would
like so long as the csv stays tidy. 
Many of the ``phippery`` API utilites,
are generalized to help you index, and otherwise
manipulate the data to you liking using any combination
of these annotations, so go wild with annotations!

.. note:: Some of the :ref:`optional workflows <sec_optional_workflows>`
    have additional required annotations, so keep an eye for those.

Keep in mind the internal datatypes are not handled perfectly
(see :ref:`a confession <sec_confession>`) 
- so it's best to keep datatypes consistant
between the columns provided. For :ref:`missing data <sec_missing_data>`, 
we reccomend empty strings, "", 
but "NaN" and "N/A" also work (hopefully) as expected.

.. todo:: reference the sample_id thing

.. _sec_input_fasta:

Sample fastq files
++++++++++++++++++

The fastq files pointed to by the sample table described above
are assumed to have uniform (trimmed) read lengths.
Note that reads are trimmed

.. todo:: Finish description
.. todo:: Should we remove the peptide length thing? is that confusing?

.. _sec_pep_anno:

Peptide table
+++++++++++++

A CSV where one of the columns must be "oligo" which
contains the oligonucleotide sequence encoding a peptide in
the phage library. 

.. todo:: Add example and explain

.. _sec_pipeline_outputs:
  
================
Pipeline results
================

The pipeline will output all results to the relative specified by the
``--dataset_prefix`` parameter.
this includes a phip_data/ directory with the pickled xarray binary file,
and optionally the tall_data/ and wide_data/ directories if specified.

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
selectio of the features found in the
phippery python API as optional during pipeline
execution. 

.. todo:: finish description

Optional Workflow: CPM Enrichment
+++++++++++++++++++++++++++++++++

``--run_cpm_enr_workflow``

.. todo:: add link to autodoc function - forgot how to do that

- help: Flag for running the enrichment workflow using counts
    per million as a pre-processing step to fold enrichment.
- wb_type: bool
- default: False

Negative Binomomial
+++++++++++++++++++

.. todo:: link to description

``--run_neg_binom_fit_predict``

- help: Flag for running negative binomial modeling
- wb_type: bool
- default: False

Z-Score
+++++++

.. todo:: link to description

``--run_zscore_fit_predict``

- help: Flag for running Z-score enrichment analysis
- wb_type: bool
- default: False

.. todo:: show example of running all the optional workflows minus virscan
    i.e. what do the sample and peptide tables look like and how long does
    take? You could poentially add some of the nextflow stats from the
    nextflow official report

VirScan Public Epitopes
+++++++++++++++++++++++

.. todo:: d

``--summarize_by_organism``

- help: Flag used to control the summary of results by organism
- wb_type: bool
- default: False

``--peptide_org_col``

- help: Column in the peptide table indicating the organism for each peptide
- wb_type: string
- default: Strain

``--peptide_prot_col``

- help: Column in the peptide table indicating the protein for each peptide
- wb_type: string
- default: Protein

``--peptide_pos_col``

- help: Column in the peptide table indicating the position within the protein for each peptide
- wb_type: string
- default: Prot_Start

``--peptide_seq_col``

- help: Column in the peptide table containing the peptide sequence (used to match against public epitopes)
- wb_type: string
- default: Prot

``--max_overlap``

- help: Maximum allowed overlap between detected peptides
- wb_type: integer
- default: 7

``--zscore_threshold``

- help: Minimum z-score threshold
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

