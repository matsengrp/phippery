.. _sec_introduction:

============
Introduction
============

.. note:: This documentation is incomplete and under development. If
    you would like to help, please open an issue or pull request at
    `GitHub <https://github.com/popgensims/stdpopsim>`_.

This is the documentation for ``phippery``, a set of functions designed to help query PhIP-Seq data.

We designed ``phippery`` to address a lack of standardized tools for analyzing counts the
*enrichment* tables obtained as a result of the protocol's sequence alignment heuristic.

First steps
-----------

 - Head to the :ref:`Installation <sec_installation>` page to get ``stdpopsim`` installed
   on your computer.


The entire pipeline
-------------------

Pipeline dependencies:  Nextflow and Docker (or have them available)
Dataset query and visualization dependencies: phippery and streamlit

we suggest:
::
  $ python -m venv phip-ex
  $ source phip-ex/bin/activate

PHIP-FLOW
^^^^^^^^^

run the nextflow pipeline locally using docker containers and native cpu
on an `example dataset <>`_ provided by the `Overbaugh Lab <>`_ from `X <publication>`_
::
  $ nextflow run matsengrp/phip-flow/PhIP-Flow.nf -C foo-phip.config -o bar.phip

PHIPPERY
^^^^^^^^

Calculate the counts per million on the raw counts
::
  $ phippery summarize bar.phip

.. note:: For a list of other 
    available commands available by the click CLI:  `$ phippery -h`  

PHIP-VIZ
^^^^^^^^
.. https://www.section.io/engineering-education/how-to-deploy-streamlit-app-with-docker/

Finally, run the phip-viz app, using the normalized sample data:
::
  $ docker pull quay.io/matsengrp/phip-viz
  $ docker run -it -v bar.phip:/ quay.io/matsengrp/phip-viz:latest
  



