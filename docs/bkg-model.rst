

.. _sec_background_modeling:

===================
Background Modeling
===================

Introduction
------------

``phippery`` provides two modeling options for estimating significance of peptide enrichment:
(1) A Gamma Poisson model as presented in `Larmen et. al. 2013 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3677742/>`_,
and (2) a Z-score method that was used in `Mina et al. 2019 <https://www.science.org/doi/10.1126/science.aay6485>`_ [#MinaMeasles]_ (and described in detail in their
`supplementary document <https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aay6485&file=aay6485_mina_sm.pdf>`_).


Z-score Method
--------------

We offer the ability to model enrichment using a Z-score method that was used in `Mina et al. 2019 <https://www.science.org/doi/10.1126/science.aay6485>`_ [#MinaMeasles]_ (and described in detail in their
`supplementary document <https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aay6485&file=aay6485_mina_sm.pdf>`_). The method uses the mock-IP
samples to bin together peptide species of similar abundance under the beads-only condition. Here, abundance can be represented in any form of normalized counts and
CPM is the default in ``phippery``. Note that the mock-IP samples are used only to determine binning.

To compute the Z-score for a peptide species in an empirical sample, identify the bin it belongs to and compute the mean and standard deviation CPM among the peptide
species in that bin. To reduce the influence of outliers, such as signal from epitope-specific binding, the highest 5% and lowest 5% of CPM values are discarded when
computing the mean :math:`\mu` and standard deviation :math:`\sigma`. Formally, for a peptide species, :math:`p`, with CPM value, :math:`n_p`, belonging to bin :math:`i`,
the Z-score is:

.. math::
	Z_p = \frac{n_p - \mu_i}{\sigma_i}



References

.. [#MinaMeasles] Mina, M.J., et al. `Measles virus infection diminishes preexisting antibodies that offer protection from other pathogens <https://www.science.org/doi/10.1126/science.aay6485>`_.
                  Science, 2019. **366** (6465): p. 599-606.
