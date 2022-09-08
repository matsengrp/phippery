
.. _sec_install_intro:

============
Installation
============


^^^^^^^^^^
Python API
^^^^^^^^^^

Installation of the ``phippery`` python package is available through

.. code-block::

  » pip install phippery

For activate development, and documentation, we recommend using the following
instructions. 

.. code-block::

  » git clone https://github.com/matsengrp/phippery.git
  » python -m venv phippery_dev_env
  » source phippery_dev_env/bin/activate
  » (cd phippery && pip install -e ".[dev]")

.. seealso:: for more information about how to contribute
  please see the :ref:`Development <sec_dev_intro>` page.


^^^^^^^^^^^^^^^^^^^
Alignments Pipeline
^^^^^^^^^^^^^^^^^^^

For running the Nextflow alignments pipeline,
we maintain a `Dockerfile`_, 
and host a pre-built `public image`_ on Quay for convinience.
with all the necessary
dependencies

To install
`docker desktop <https://www.docker.com/products/docker-desktop>`_ 
(Desktop is just fine), 

and 
`Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_. 
Be sure to have the command line tools for both.

To install `docker` on most unix OS:

.. code-block:: bash

    curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

To test the `docker` install:

.. code-block:: bash

   » docker -v
   Docker version 20.10.1, build 831ebea
   » docker run docker/whalesay cowsay boo

To install `Nextflow` on most Unix OS

.. code-block:: bash

   » curl -s https://get.nextflow.io | bash 

To test the `Nextflow` install:

.. code-block:: bash

   » nextflow -v
   nextflow version 20.04.1.5335

.. note:: The pipeline has been tested up to Nextflow version 22.04.3.5703   

.. _sec_installation_phippery:

^^^^^^^^^^^^^^^^^^^^^^^^
Python CLI/API
^^^^^^^^^^^^^^^^^^^^^^^^

The API and CLI for dataset querying and for visualization 
are implemented in python, with commonly used data science
dependencies such as 
`numpy <https://numpy.org/doc/stable/user/basics.dispatch.html>`_ and
`pandas <https://pandas.pydata.org/>`_, 
along with a few other less traditional dependencies such as 
`xarray <http://xarray.pydata.org/en/stable/>`_ and
`streamlit <https://docs.streamlit.io/en/stable/>`_. 
For a full list of respective dependencies, see the 
`phippery <https://github.com/matsengrp/phippery/blob/master/requirements.txt>`_ and 
`phip-viz <https://github.com/matsengrp/phip-viz/blob/main/requirements.txt>`_ 
requirements files. We maintain and provide
`pip <https://pypi.org/>`_ 
installation options for either of these tools, or a 
`Docker <https://www.docker.com/>`_ 
image with dependencies 
(useful for extending the PhIP-Flow pipeline).


pip
^^^

Currently you need to do `clone` and install of the dependencies locally,
keeping each of the directories separate.
 
.. note::
   we suggest using pip + venv to install ``phippery`` and ``phipviz``

   .. code-block::

     »  python -m venv phip-env
     » source phip-env/bin/activate

.. code-block::     

   » git clone https://github.com/matsengrp/phippery.git
   » (cd phippery && pip install .)

^^^^^^^^^^^^^
Streamlit app
^^^^^^^^^^^^^

.. code-block::  

  » git clone https://github.com/matsengrp/phip-viz.git
  » (cd phip-viz && pip install -r requirements.txt)

.. note:: phippery PyPI NOT RELEASED, YET. Coming soon



