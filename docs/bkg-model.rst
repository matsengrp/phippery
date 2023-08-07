

.. _sec_background_modeling:

===================
Background Modeling
===================

Introduction
------------

In addition to running the alignments and computing the counts, ``phippery`` also provides workflows 
to model the background with mock-IP samples and compute non-specific peptide enrichments.
By default the ``phip-flow`` pipeline will run the
`edgeR <https://doi.org/10.1093%2Fbioinformatics%2Fbtp616>`_ [#edgeR2010]_ workflow as described in the 
`Chen et al. 2022 <https://doi.org/10.1093/bioinformatics/btac555>`_ [#ChenBEER]_ paper.
Optionally, we also provide a simpler Z-score method to evaluate the significance of peptide enrichment relative to
background that was used in `Mina et al. 2019 <https://www.science.org/doi/10.1126/science.aay6485>`_ [#MinaMeasles]_. Each is described in greater detail below. Note that these workflows are not mutually exclusive, (i.e. do not overwrite the counts, cpm, or any other default/optional outputs from the pipeline. You may run none or all of the workflows in tandem, then do with the :ref:`combined results <sec_pipeline_outputs>` as you wish.

.. _sec_edger:

edgeR/BEER Method
-----------------
`Chen et al. 2022 <https://doi.org/10.1093/bioinformatics/btac555>`_ adapts the `edgeR <https://doi.org/10.1093%2Fbioinformatics%2Fbtp616>`_ tool to compute
fold-change with respect to mock-IP samples and p-values of peptide enrichment. Optionally, you may run the 
`BEER (Bayesian Estimation Enrichment in R) method <https://bioconductor.org/packages/release/bioc/vignettes/beer/inst/doc/beer.html#beer-bayesian-estimation-enrichment-in-r>`_,
which is statistically more powerful and may be better at identifying significantly enriched peptides with lower fold-changes. 
The trade-off for using the BEER method is longer run-time.
By default, the ``phip-flow`` pipeline runs EdgeR, but not BEER. 
see :ref:`Optional Parameters in the pipeline documentation <sec_optional_workflows>` for more. 

Z-score Method (optional)
-------------------------

``phippery`` can also optionally run the Z-score method to compute the significance of peptide enrichment relative to background.
This Z-score method used in `Mina et al. 2019 <https://www.science.org/doi/10.1126/science.aay6485>`_ is described in detail in their
`supplementary document <https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aay6485&file=aay6485_mina_sm.pdf>`_. The method takes the mock-IP
samples to bin together peptide species of similar abundance under the beads-only condition. Here, abundance can be represented in any form of normalized counts and
CPM is the default in ``phippery``. Note that the mock-IP samples are used only to determine binning.

To compute the Z-score for a peptide species in an empirical sample, identify the bin it belongs to and compute the mean and standard deviation CPM among the peptide
species in that bin. To reduce the influence of outliers, such as signal from epitope-specific binding, the highest 5% and lowest 5% of CPM values are discarded when
computing the mean :math:`\mu` and standard deviation :math:`\sigma`. Formally, for a peptide species, :math:`p`, with CPM value, :math:`n_p`, belonging to bin :math:`i`,
the Z-score is:

.. math::
	Z_p = \frac{n_p - \mu_i}{\sigma_i}

References
----------

.. [#edgeR2010] Robinson, M.D., McCarthy, D.J., and Smyth, G.K.
                `edgeR: a Bioconductor package for differential expression analysis of digital gene expression data <https://doi.org/10.1093%2Fbioinformatics%2Fbtp616>`_.
                Bioinformatics, 2010. **26** (1): p. 139-140.

.. [#ChenBEER] Chen, A., et al. `Detecting and quantifying antibody reactivity in PhIP-Seq data with BEER <https://doi.org/10.1093/bioinformatics/btac555>`_.
               Bioinformatics, **38** (19): p. 4647-4649.

.. [#MinaMeasles] Mina, M.J., et al. `Measles virus infection diminishes preexisting antibodies that offer protection from other pathogens <https://www.science.org/doi/10.1126/science.aay6485>`_.
                  Science, 2019. **366** (6465): p. 599-606.
