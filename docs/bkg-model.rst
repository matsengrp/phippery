

.. _sec_background_modeling:

===================
Background Modeling
===================

Introduction
------------

Ideally, the phages that are pulled down in the immunoprecipitate are bound to antibodies, but in practice some unbound phages are
also attracted by the magnetic beads. This introduces a background of peptide counts. To assess this background, PhIP-Seq runs are
performed without the presence of antibodies; these *mock-IP* samples are used to derive a background model for each peptide 
species in the library. Peptide counts from empirical samples can then be compared with the associated background model to compute 
a p-value -- how likely the observed count is a background interaction. Therefore, genuine enrichment from antibody-peptide binding 
should result in a small p-value. In ``phippery``, the negative binomial model is implemented as a fit per peptide species over 
mock-IP data.


Normalization via Size Factors
------------------------------

Before fitting to the mock-IP data, counts must be normalized to account for variations in total reads across samples.
It might be tempting to turn to the commonly used *counts per million* (CPM) normalization, but this is problematic if the
total reads of the samples are very different from one million. This is because read counts are discrete, not continuous,
and scaling the data arbitrarily can yield a distribution that cannot be simply described by a probability model. 
To choose a normalization factor that is more natural to the dataset, we employ the *size factor* quantity [#SizeFactors]_:

- A *pseudo-reference sample* is derived from computing the geomtric mean across all samples for each peptide
- For a specific sample to be normalized, consider for each peptide the ratio of the sample count to the pseudo-reference count; the median of these ratios is the normalization factor.

This method is designed to be less influenced by outlier counts from 
very highly enriched peptides. Both mock-IP and empirical samples should be considered together when
normalizing via size factors. In ``phippery``, the ``size_factor()`` function generates a new data table of size factor
normalized counts based on the samples in the input ``xarray.Dataset``.

As an example, consider the distribution of counts for a peptide across 207 mock-IP samples,

.. image:: images/counts.png
	:width: 500
	:align: center
	
The average total reads in these mock-IP samples is about 200 000. Consequently, the CPM distribution contains many "gaps"
and is difficult to model. (The range of the histogram below is limited to 200 to show the problematic gaps, but the full
distribution extends much further.)

.. image:: images/cpm_norm.png
	:width: 500
	:align: center

With size factor normalized counts, the distribution is more reasonable to model.

.. image:: images/sf_norm.png
	:width: 500
	:align: center


Negative Binomial Model
-----------------------

Peptide counts generally show over-dispersion with respect to the ideal Poisson distribution, and the negative binomial
distribution is a commonly used model to address this situation. A model is fit to the mock-IP data for each peptide species
and two parameters are returned, corresponding to the 
`SciPy <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html>`_ 
parametrization of the negative binomial. Occasionally, a fit may fail to converge, usually due to little or no representation
in the mock-IP samples. The table below lists the possible outcomes of a fit.

.. list-table:: Return values for negative binomial fit
   :widths: 50 50
   :header-rows: 1

   * - Returns
     - Reason
   * - :math:`\mbox{size}\geq0, \mbox{prob}\geq0`
     - Fit converged normally
   * - :math:`\mbox{size}=\mbox{prob}=-1`
     - All counts are zero
   * - :math:`\mbox{size}=\mbox{prob}=-2`
     - Fit failed to converge

Below are some examples of fits to data.

.. image:: images/SUB2_nb_fit_peptide_2069.png
	:width: 500
	:align: center
.. image:: images/MEGSUB_nb_fit_peptide_6149.png
	:width: 500
	:align: center
.. image:: images/MEGSUB_nb_fit_peptide_10035.png
	:width: 500
	:align: center


Z-score Method
--------------

If there are few mock-IP samples (~10 or less) available, the negative binomial model fits may struggle to converge. An alternative method implemented in ``phippery``
is a Z-score method that was used in `Mina et al. 2019 <https://www.science.org/doi/10.1126/science.aay6485>`_ [#MinaMeasles]_ (and described in detail in their
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


.. [#SizeFactors] Anders, S. and Huber, W., `Differential expression analysis for sequence count data
                  <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106>`_. Genome Biology, 2010. **11**:R106.

.. [#MinaMeasles] Mina, M.J., et al. `Measles virus infection diminishes preexisting antibodies that offer protection from other pathogens <https://www.science.org/doi/10.1126/science.aay6485>`_.
                  Science, 2019. **366** (6465): p. 599-606.