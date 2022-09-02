
.. _sec_streamlit_app:

=================
Streamlit Viz App
=================

For convenience, we have created a simple 
`streamlit <https://streamlit.io/>`_ application for viewing a ``phippery``
dataset (i.e. "data.phip") as 
`altair heatmaps <https://altair-viz.github.io/gallery/simple_heatmap.html>`_.
The application allows you to subset you data based upon annotations that 
were in the sample or peptide tables (see :ref:`pipeline inputs <sec_pipeline_inputs>`).

To use the application on the :ref:`example Pan-CoV dataset <sec_clone_template>`,
follow the :ref:`installation instructions <sec_install_intro>`.
Next, move into the directory where ``phip-viz`` was cloned, and link the
dataset to be in the same directory.

.. code-block:: bash

    $ cd phip-flow
    $ ln -s ../phip-flow/results/pickle_data/data.phip ./

.. note:: This command assumes you cloned both the phip-viz, phippery, and phip-flow
    repositories into the same parent directory.

Next, use the ``streamlit run`` command to launch the application.

.. code-block::

    $ streamlit run streamlit_app.py

Your default browser will then open the application and look something like this

.. figure:: images/launch_viz_app.png
  :width: 700
  :alt: example results
  :align: left

As a quick example, click on the expander drop bar on the right side, labeled 
"Working Peptides" under ``Peptide table``,
This will show us a brief overview of the peptides in our dataset for reference, or to 
add queries that will narrow the scope of our visualizations to a smaller subset.
Notice that the application gives some helpful suggestions based upon the data available.

We'll start by looking at only peptides from the spike protein of SARS-CoV-2 and 229E 
viruses. To do this, type ``(Virus in ['SARS', '229E']) & (Protein == 'spike')`` in the entry box
denoted ``Peptide Query Condition``, and hit the ``<Enter>`` key. You'll 
notice the query and it's unique key (denoted starting with "q<n>") 
was added to the table. The query was applied to the working
dataset and you can now see the summary of peptides has changes accordingly.

.. figure:: images/query_peptides_viz_app.png
  :width: 700
  :alt: example results
  :align: left

Next, well scroll down to the ``Visualize Enrichment Heatmap`` section to look at our data.
Enter the following options to view the enrichments aggregated by sample type and locus.

.. figure:: images/heatmap_viz_app.png
  :width: 700
  :alt: example results
  :align: left

And that's it! For other features such as saving images, uploading query tables, and more, 
click on and of the ``?`` boxes to explore the various options.
