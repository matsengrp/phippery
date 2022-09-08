
.. _sec_install_intro:

============
Installation
============

.. _sec_installation_phippery:

^^^^^^^^^^
Python API
^^^^^^^^^^

Installation of the ``phippery`` python package is available through

.. code-block::

  » pip install phippery

.. seealso:: for more information about how to contribute
  please see the :ref:`Development <sec_dev_intro>` page.


^^^^^^^^^^^^^^^^^^^
Alignments Pipeline
^^^^^^^^^^^^^^^^^^^

For running the ``Nextflow`` alignments pipeline,
we maintain a 
`Dockerfile <https://github.com/matsengrp/phip-flow/blob/main/docker/Dockerfile>`_, 
and host a pre-built 
`public image <https://quay.io/repository/jgallowa/phip-flow>`_ 
on quay.io for convenience.

.. note::
    Note that the quay image is pointed to by default
    in the pipeline configuration file --
    meaning if you have both ``Docker`` and ``Nextflow``
    installed, you have everything you need to run the pipeline 
    with ``--profile docker`` specified after the ``nextflow run``
    command. The other option would be to pull, 
    or build the image yourself, then run the pipeline
    within that container.

To install
`docker desktop <https://www.docker.com/products/docker-desktop>`_ 
(Desktop is just fine), 
on most unix OS:

.. code-block:: bash

   » curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

To test the `docker` install:

.. code-block:: bash

   » docker -v
   Docker version 20.10.1, build 831ebea
   » docker run docker/whalesay cowsay boo

To install
`Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_
on most Unix OS:

.. code-block:: bash

   » curl -s https://get.nextflow.io | bash 

To test the `Nextflow` install:

.. code-block:: bash

   » nextflow -v
   nextflow version 20.04.1.5335


^^^^^^^^^^^^^
Streamlit App
^^^^^^^^^^^^^

.. code-block::  

  » git clone https://github.com/matsengrp/phip-viz.git
  » (cd phip-viz && pip install -r requirements.txt)
