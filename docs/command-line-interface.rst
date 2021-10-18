.. _sec_cli_intro:

======================
Command Line Interface
======================

For all types of analysis outside of read alignment and visualization, 
we recommend using the Command Line Interface (CLI) accessed using the `phippery` command.
The CLI is written using the 
`Click <https://click.palletsprojects.com/en/8.0.x/>`_
library, and thus both `phippery -h`, and `phippery COMMAND -h` will provide
the same information provided below

.. note:: For specific usages of these commands please see TODO=

.. warning:: This CLI is still under activate development. 
  we appreciate your feedback and patience

.. click:: cli:cli
   :prog: phippery
   :nested: full

.. click:: cli:phipflowcli
   :prog: phipflow
   :nested: full
