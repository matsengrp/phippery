.. _sec_introduction:

============
Introduction
============

.. note:: This documentation is incomplete and under development. If
  you would like to help, please open an issue or pull request at
  `GitHub <https://github.com/matsengrp/phippery/>`_.

Welcome to the documentation for ``phippery``. 
Here, we present a set of data analysis tools designed to aid in the
exploration of `antibody binding <TODO>`_ studies.
We have designed a suite software to be flexible to the data resulting from both
*common* and *highly customized* phage-display libraries commonly used in
`Phage Immuno-Precipitation (PhIP-Seq) <https://www.nature.com/articles/s41596-018-0025-6>`_ experiments.

===============================================

.. figure:: images/phippery-suite-5.svg
  :width: 1000
  :alt: Alternative text
  :align: center

  **software tools schematic:** A cartoon flow 
  chart giving a description of the workflow using
  each of the the tools presented here. (A) Shows
  the input from the researcher; First, a set of 
  demultiplexed files for the samples of interest, 
  as well as a sample annotation table describing
  features relevent to a study, as well as the file
  paths to each of the respective fastq files.
  Additionally, there is a peptide table which
  also describes the peptide features of interest
  (virus, protein, locus, etc), and importantly
  the reference oligo sequence the pipeline should
  align sample reads to. (B) Visualizes the set of
  CLI tools for data formatting, storing, transforming, 
  fitting models to get estimates of significance.
  Most operations result in "layering" specific
  transformations onto the primary data structure.
  For more detail on code structure, and the 
  available command line interface (CLI), see the
  :ref:`under the hood <sec_python_intro>` section.
  Finally, (C) takes the layered data structure
  and allows users to aggregate and visualize 
  any combination of the sample or peptide
  annotation features provided. The final export
  can be either the visualization itself, or the
  underlying raw data (tall or wide) for plotting 
  with any of your favorite packages.

===============================================

+++++++++++++++
Getting Started
+++++++++++++++


- Head over to the :ref:`Installation <sec_install_intro>` 
  page to get installation instructions for each of the three tools described above.

- To get a feel for running each of the three related tools pictured above, 
  we suggest following a 
  :ref:`walk through <sec_quick_start>` of running the alignments, CLI, and visualization
  app on some empirical data from `Stoddard et al. 2021 
  <https://www.cell.com/cell-reports/fulltext/S2211-1247(21)00506-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124721005064%3Fshowall%3Dtrue>`_: 
  *Epitope profiling reveals binding signatures of SARS-CoV-2 immune response 
  in natural infection and cross-reactivity with endemic human CoVs* 

- If the Take a look at the :ref:`CLI commands <sec_cli_intro>` offered by the `phippery` command 
  to perform you own analysis steps. 
  

- For creating your own dataset annotations, configuration, and NGS file structure
  first, Review the 
  :ref:`Annotations <sec_pipeline_anno>` section, 
  as many of the tools here require 
  these tables as input for analyzing your own data.

++++++++++
Background
++++++++++

Since the formation of PhIP-Seq, `many <TODO>`_
studies have detailed important insight into the antibody 
(humoral) response of individuals
under various conditions usually involving an 
`antigenic <TODO>`_-provoked adaptive immune response.
In brief detail, the advent of modern 
`oligonucleotide synthesis <TODO>`_
allows researchers to generate *highly* multiplexed assays, with up to the order
:math:`10^{4}` peptides being expressed using `phage display <TODO>`_ libraries.
As a result, the PhIP-Seq protocol uses modern techniques
to quantify the concentration (enrichment), of antibody-peptide binding 
events when serum is extracted and presented to the short linear proteins
One impressive and commonly used phage libraries is the 
`Virscan <TODO>`_ library -- with <X> peptides generated from <X> unique Virus-Proteins.
While impressive, researchers have started creatively developing 
smaller, and more custom libraries to explore even more nuanced questions concerning exact
`linear epitope` intervals, as well as exploring differential selection of those epitopes
under various conditions of mutations. 

With rapid progress, it's not surprising that researchers
find themselves either piecing together published code with their own, 
or writing completely novel analysis scripts on the fly. 
While impressive, this tends to make comparison of analysis
between the many new studies quite difficult to interpret and apply to ones own data.
The goal here is to provide some *efficient* and *unit-tested*
infrastructure for; computing enrichment 
(given some demultiplexed `NGS data files <TODO>`_), 
data formatting, storing, transforming, fitting models to
PhIP-Seq data deriving from both _common_ libraries, to completely novel custom libraries.
Each of the tools presented here can be used separately, or in conjunction for the
rapid exploration of PhIP-Seq data.

+++++++++++++++++++++++++++++
Licensing and Acknowledgement
+++++++++++++++++++++++++++++

This work is provided by members of the 
`Matsen <TODO>`_ and 
`Overbaugh <TODO>` groups at the
`Fred Hutchinson Cancer Research Center <TODO>`_.
The software is publicly available and licensed under the 
`GNU GENERAL PUBLIC LICENSE <TODO>`_
The work presented is funded by the **NIH**, **NSF**, and **HHMI**.

For questions or concerns about these using tools, 
feel free to email jgallowa (at) fredhutch
If you find these tools useful for 
your own research studies, please cite <X>


