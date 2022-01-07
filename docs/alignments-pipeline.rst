
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

.. _sec_sam_anno:

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

