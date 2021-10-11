

.. _sec_escape_profile_comparisons:

=========================
Comparing Escape Profiles
=========================

Introduction
------------

Phage-DMS [1] is an application of PhIP-Seq for deep mutational scanning (DMS).
The Phage-DMS method measures changes in antibody binding under mutations to 
viral proteins. Of particular interest are mutations that cause significant loss
of binding, potentially allowing the virus to escape the immune response. When 
considering only mutations that cause single amino acid substitutions, the results 
can be summarized in an escape profile showing the impact on binding for every
possible substitution at each location.

This page describes a method for quantifying the comparison of escape profiles 
(e.g. between two individuals). There are two main steps:

1. Assigning a similarity score for each site quantifying the comparison of differential selection between the two samples.
2. Computing a weighted sum of similarity scores across the sites of interest to obtain a single score for the region.


Similarity score at a site
--------------------------

The task of comparing two escape profiles at a site is handled in the framework of 
an optimal transport problem. In brief, an optimal transport problem involves two
distributions, a specification of the cost to transport "contents" between them, 
and the objective of determining how to transform one distribution into the other
while minimizing the total cost. This minimum cost reflects the level of similarity
between the two distributions; the lower the cost, the more similar they are.

Each escape profile is represented as a binned distribution. There are 20 amino acids
but a mutation could result in negative or positive differential selection, respectively
a loss or a gain in binding, which are two different effects. Therefore two bins are
associated with each amino acid, one for negative differential selection and one for
positive differential selection, resulting in a 40-bin distribution. The bin contents
are the relative contribution to the total absolute differential selection at the site.
For any escape profile, at most 19 bins are non-zero because an amino acid substitution
cannot contribute negative and positive differential selection simultaneously, and the 
wild type amino acid contributes zero by definition.

The cost matrix is based on the BLOSUM62 substitution matrix [2], which assigns an integer
score between two amino acids based on how often the substitution is observed in empirical
data. Hence, the entries of the BLOSUM62 matrix, denoted here by :math:`M_{ij}`, quantifies
an overall similarity between amino acids :math:`i` and :math:`j`. Using the BLOSUM62 matrix
allows for escape profile comparisons beyond just contributions from matching amino acids 
in the two samples. For transport between like-sign differential selection (e.g. negative
with negative) "contents" between :math:`i` and :math:`j`, the cost is given by,

.. math::
	C_{ij} = \exp\left(-M_{ij}/7\right).

The exponentiation is applied so that cost values are positive. Higher values in the BLOSUM62
matrix corresponds to higher similarity between amino acids, and vice versa, so the negative 
sign is applied to assign the corresponding cost. The factor of 7 is included so that the 
range of cost values lie within a factor of 10, to avoid a single site from overwhelming
contributions from all other sites when we aggregate across a region.

For transport between opposite-sign differential selection, we disregard similarity between
amino acids and simply assign a fixed cost of :math:`C_{\mbox{max}} = \exp\left(4/7\right)`,
which is the maximum cost possible with the BLOSUM62 matrix.

Putting this altogether, the complete cost function is a :math:`40\times40` matrix that can
be expressed in blocks of :math:`20\times20` sub-matrices as the following: 

.. math::

	\begin{bmatrix}
		C_{ij} & C_{\mbox{max}} \\
		C_{\mbox{max}} & C_{ij}
	\end{bmatrix},

where the off-diagonal blocks are sub-matrices with :math:`C_{\mbox{max}}` for all entries.

To solve the optimal transport problem, we use the Python Optimal Transport package [3].
Because the obtained minimum cost is inversely related to how similar the two escape profiles
are, we define the *similarity score* to be the reciprocal of this cost value. Our interest
leans more towards identifying escape profiles that are consistent, and working with
similarity score rather than cost makes the interpretation a little easier when we aggregate
scores across sites in a region. Namely, we attribute a high regional similarity score to
having several sites with high similarity, as oppose to a low regional cost due to lacking
sites with high cost. Of course, the conclusions won't change either way; cost and similarity
are just two sides of the same coin.


Similarity score for a region
-----------------------------

Having defined and computed the similarity score for each site, a score for the comparison
across a region of sites is calculated by a weighted sum of similarity scores. The weights
are assigned so that sites contributing greatly to the overall escape in the region in
both profiles are given the most importance.

First, compute the relative contribution of each site to each profile. This is the summed
absolute scaled differential selection at a site divided by the sum over all sites in the
region. Denote these relative contributions by :math:`\left\{\alpha_k\right\}` and 
:math:`\left\{\beta_k\right\}` for the two profiles, where :math:`k` runs over the sites.
Then for each site, we choose the minimum,

.. math::

	\mu_k = \min\left(\alpha_k,\,\beta_k\right).

The motivation to take the minimum (as oppose to the average) is because we wish to lower
the importance of the similarity score at a site that has very different relative contributions
in the two profiles (i.e. very high in one and very low in the other).

Finally, we normalized so that the sum of weights is 1,

.. math::
	
	w_k = \frac{\mu_k}{\sum_k \mu_k}.

Given the set of similarity scores across the sites, :math:`\left\{s_k\right\}`, the similarity
for the region is,

.. math::

	S = \sum_k w_k s_k


References
----------

[1] Garrett, M.E., et al., `Phage-DMS: A Comprehensive Method for Fine Mapping of Antibody Epitopes <https://doi.org/10.1016/j.isci.2020.101622>`_. iScience, 2020. 23(10): p. 101622.

[2] Henikoff, S., Henikoff, J.G., `Amino Acid Substitution Matrices from Protein Blocks <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/>`_. PNAS, 1992. 89(22): p. 10915-10919.

[3] Flamary, R., et al., `POT: Python Optimal Transport <https://pythonot.github.io/>`_. JMLR, 2021. 22(78): p. 1-8.