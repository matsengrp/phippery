
.. _sec_dev_intro:

===========
Development
===========

Developer Install
^^^^^^^^^^^^^^^^^

For activate development, and documentation, we recommend using the following
instructions. 

::

  » git clone https://github.com/matsengrp/phippery.git
  » python -m venv phippery_dev_env
  » source phippery_dev_env/bin/activate
  » (cd phippery && pip install -e ".[dev]")

Next, run the tests to make sure everything is working properly.

::

  » pytest -vv


Building Documenation
^^^^^^^^^^^^^^^^^^^^^

To edit the documentation seen here,
simply edit the respective `.rst` file 
(following the git workflow described below) 
in the `docs/` subdirectory. Once edited, you can check 
the edits are rendered correctly by building the docs locally

.. code-block::

  » cd docs/
  » make clean && make html

Then open the index file build at `_build/html/index.html`
with a browser of choice to inspect changes.

Once the changes have been approved and merged into the main branch
the documentation will automatically build and deploy.

  
Contributing
^^^^^^^^^^^^

TODO

