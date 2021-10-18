
.. _sec_install_intro:

============
Installation
============

Each of the tools presented here can be run with the installation of
`Docker`, `Nextflow`, `pip`, and `git`. 
The details of installation and updating each are described
below. For a soup-to-nuts example of running all three tools together, see
:ref:`Quick Start <sec_quick_start>` section.

.. note:: The software presented here is still under construction and 
    considered to be in the "Beta" stage of production. 
    Please expect and check back for innevitable changes, 
    for questions and/or suggestions, please feel welcome 
    to contact jgallowa (at) fredhutch.org

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Alignment pipeline dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, make sure to have a working install of
`docker desktop <https://www.docker.com/products/docker-desktop>`_ 
(Desktop is just fine) and 
`Nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_. 
Be sure to have the command line tools for both.

To test the `docker` install:

.. code-block:: bash

   » docker -v
   Docker version 20.10.1, build 831ebea
   »  docker run docker/whalesay cowsay boo

To test the `nextflow` install:

.. code-block:: bash

   » nextflow -v
   nextflow version 20.04.1.5335
   » nextflow run hello

From here, we can simply use `Nextflow's git aware <TODO>`_ 
infrastructure to run the bleeding edge script directly from the source 
`git repository <https://github.com/matsengrp/phip-flow>`_.
For example, we are now ready to run the pipeline like so,
given a config file ``foo.config``

.. code-block:: bash

   » nextflow -C foo.config run matsengrp/phip-flow/PhIP-Flow.nf

Of course, this is assuming you've got all the configuration
files ready. For a quick introduction to the input files
with some examples, check out the :ref:`Examples <sec_quick_start>`
page. For even more details on input formatting and preparing
to create and run your own pipeline, please see the
:ref:`Alignments Pipeline <sec_pipeline_intro>` page.

.. tip:: If you would like to retain a copy of the Nextflow 
  script locally for modification, simply clone 
  the source code `pipeline repository <TODO>`_, 
  or use the `--recurse-submodules` flag when cloning 
  the template as described in the 
  :ref:`Example run <sec_clone_template>` section of the
  full Pan-CoV example pipeline walk through.

.. _sec_installation_phippery:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
phippery \& phip-viz dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The command line interfaces for dataset querying and for visualization 
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

we suggest using pip + venv to install ``phippery`` and ``phipviz``
::

  $ python -m venv phip-ex
  $ source phip-ex/bin/activate
  $ pip install phippery phipviz


.. note:: PyPI NOT RELEASED, YET.
  Currently you need to do this:
  ``pip install git+https://github.com/matsengrp/phippery.git@52b8c5fcd0c4c727fe760b17a3820a60eada2bf3``


Docker
^^^^^^

We also provide a full container image with everything you need to
run both ``phippery`` and ``phipviz`` 

.. note:: Docker image NOT RELEASED, YET.

Developer Install
^^^^^^^^^^^^^^^^^

For activate development, and documentation, we recommend using the following
instructions. 

::

  git clone https://github.com/matsengrp/phippery.git
  python -m venv phippery_dev_env
  source phippery_dev_env/bin/activate
  pip install -e ".[dev]"

.. seealso:: for more information about how to contribute
  please see the :ref:`Development <sec_dev_intro>` page.
