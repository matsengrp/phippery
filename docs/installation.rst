.. _sec_installation:

============
Installation
============

Each of the tools presented here can be run with the installation of
``Docker``, ``Nextflow``, ``pip``, and ``git``. The details for each are described
below. For a soup-to-nuts example of running all three tools together, see
:ref:`Quick Start <sec_quick_start>` section.



.. _sec_installation_phip_flow:

phip-flow
+++++++++

If interested in using our ``Nextflow`` alignment pipeline to process your demultiplexed
sample  files 

To install ``Nextflow``, see their 
`instructions <https://www.nextflow.io/docs/latest/getstarted.html#installation>`_.
Once installed, clone the 
`phip-flow repository <https://github.com/matsengrp/phip-flow>`_. 

.. note:: If you have a config file already specifying 
  with and your own data, you can run the bleeding edge
  pipeline script using Nextflow's 
  `git aware infrastructure 
  <https://www.nextflow.io/docs/latest/sharing.html#running-a-pipeline>`_:
  ::
    
    nextflow run matsengrp/phip-flow/PhIP-Flow.nf TODO - Finish

.. _sec_installation_phippery:

phippery \& phip-viz
++++++++++++++++++++

Both the command line interface for queries on the dataset, and visualization, are implemented 
in python, mostly depending on common data science dependencies such as 
`numpy <https://numpy.org/doc/stable/user/basics.dispatch.html>`_ and
`pandas <https://pandas.pydata.org/>`_, 
along with a few other less traditional dependencies such as 
`xarray <http://xarray.pydata.org/en/stable/>`_ and
`streamlit <https://docs.streamlit.io/en/stable/>`_. 
For a full list of respective dependencies, see the 
`phippery <https://github.com/matsengrp/phippery/blob/master/requirements.txt>`_ and 
`phip-viz <https://github.com/matsengrp/phip-viz/blob/main/requirements.txt>`_ requirements files.
We maintain and provide
`pip <https://pypi.org/>`_ 
installation options for either of these tools, or a 
`docker <https://www.docker.com/>`_ 
image with dependencies for both.


pip
^^^

we suggest using pip + venv to install ``phippery`` and ``phipviz``
::

  $ python -m venv phip-ex
  $ source phip-ex/bin/activate
  $ pip install phippery phipviz


.. note:: PyPI NOT RELEASED, YET.
  Currently you need to do this:
  `pip install git+https://github.com/matsengrp/phippery.git@52b8c5fcd0c4c727fe760b17a3820a60eada2bf3`


docker
^^^^^^

We also provide a full container image containing everything you need to
run both ``phippery`` and ``phipviz`` 

.. note:: Docker image NOT RELEASED, YET.

