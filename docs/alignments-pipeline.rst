
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

The ``PhIP-Flow`` pipeline :ref:`inputs <sec_pipe_inputs>` are 
just two CSV files with only a single column
requirement in each. 
Concretely, the pipeline requires that a user specifies; 
(1) a sample annotation table, which (at a minimum) must include a column header"*fastq_filename*",
with each respective row providing a path (relative to workflow launch) to each sample in the study.
(2) a peptide annotation table, which (at a minimum) must include the column header, "*oligo*",
with each of it's respective rows providing the oligonucleotide encoding for each peptide in the
phage library used.
The default workflow then performs all of the major steps in processing the raw data and 
obtaining a enrichment dataset (along with some other statistical goodies).
The pipeline will output a pickle dump'd ``Xarray DataSet``, or optionally
two common CSV formats, tall and wide for the user to query with 
their own favorite analysis tools.

To be concise, the processing steps were as follows;
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
(5) It then compute some basic stats to
include along with the raw alignment counts
before it


Quickstart 
^^^^^^^^^^

To install `docker` on most unix OS:

::

    $ curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

To install `Nextflow` on most Unix OS

::

    $ curl -s https://get.nextflow.io | bash 

Run the example Pan-CoV data

::

    $ nextflow run matsengrp/phip-flow -profile docker

^^^^^^^^^^^
Input files
^^^^^^^^^^^
.. _sec_sam_anno:

Sample Table 
++++++++++++

TODO

.. _sec_pep_anno:

Peptide Table
+++++++++++++

TODO

^^^^^^^^^^^^^^^^^^^
Pipeline parameters
^^^^^^^^^^^^^^^^^^^

TODO

`--results` 
+++++++++++
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to CalliNGS-NF's folder: `results` 

Example: 
::

    $ nextflow run matsengrp/phip-flow --results /home/user/my_results

  
^^^^^^^^^^^^^^^^
Pipeline results
^^^^^^^^^^^^^^^^

TODO

^^^^^^^^^^^^
Requirements 
^^^^^^^^^^^^

TODO

^^^^^^^^^^^^^^^^^^^^^^
Directed Acyclic Graph
^^^^^^^^^^^^^^^^^^^^^^

.. image:: images/dag.svg
  :width: 600
  :alt: Alternative text
  :align: center
 
^^^^^^^^^^
Components
^^^^^^^^^^

TODO

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Licensing and Acknowledgement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Huge thanks to all the folks at 
`CalliNGS-NF <https://github.com/CRG-CNAG/CalliNGS-NF>`_ 
For inspiration while structuring this pipeline.

This work is provided by members of the 
`Matsen <https://matsen.fredhutch.org/>`_ and 
`Overbaugh <https://research.fredhutch.org/overbaugh/en.html>`_ groups at the
`Fred Hutchinson Cancer Research Center <https://www.fredhutch.org/en.html>`_.
The software is publically available licenced under the 
`GNU GENERAL PUBLIC LICENSE <https://opensource.org/licenses/gpl-license.php>`_.
The work presented is funded by the **NIH**, **NSF**, and **HHMI**.

For questions or concerns about these using tools,
feel free to email jgallowa (at) fredhutch
If you find these tools useful for your own research studies, please cite <X>

