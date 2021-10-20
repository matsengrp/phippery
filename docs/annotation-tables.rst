
.. _sec_anno_intro:

=================
Annotation Tables
=================

.. _sec_pipeline_anno:

Informative annotations for the enrichment matrix rows (peptides), 
and columns (samples) are at the core of our approach to analysis of this
type of data. The high-throughput, brute-force nature of this protocol
results in what may be most understand as *many* small experiments 
- providing us with information about *binding profiles* for each 
of serum samples in a study. 
These binding profiles are comprised of tens, to hundreds-of-thousands
of measurements, each themselves representing a potential binding event to 
one of the peptides in the entire phage display library.

Because this information is so useful for computing various transformations on the data,
and so useful for keeping large datasets organized and coherent,
we require the user provides these annotations tables in a strait-forward, 
albeit *specific* format. 
If you run the 
:ref:`alignment pipeline <sec_pipeline_intro>`
provided here, then this format is prepared and currently pickle 
dumped to a binary file, usually with a ".phip" extension.

While there are really no
*required* annotations for this format,
but the tools described here are pretty 
useless without data describing groups effectively in the dataset.
We describe 
:ref:`common feature groups <sec_sam_anno>`:
below - but note that many
functions have generalized capability
to facilitate creative analysis tailored
to specific study.


.. figure:: images/xarray-format.svg
  :width: 400
  :alt: phippery xarray format
  :align: left

  **Xarray Dataset Format** A cartoon representation
  of the format output by the pipeline for 
  we requre for using phippery functions.
  Concretely, for a matrix, `\mathcal{M}_{i}{j}` 


.. _sec_sam_anno:

++++++++++++
Sample Table 
++++++++++++


+++++++++++++
Peptide Table
+++++++++++++



++++++++++
Discussion
++++++++++

Missing Data
------------


When dealing with missing values in the annotation tables, we use the 
`pd.convert_dtypes <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.convert_dtypes.html>`_
function to best allow for missing annotations, while maintaining the integrity of
the inferred datatype. It is highly recommended you stay consistent with datatype for feature annotations,
i.e. try not to mix values like `1` (integer), `6.7` (float), and `hello_world` (string) in any one of the columns. 
For missing data of any type, 
The following values will be interpreted as `NaN`:

::

  ‘’, ‘#N/A’, ‘#N/A N/A’, 
  ‘#NA’, ‘-1.#IND’, ‘-1.#QNAN’, ‘-NaN’, ‘-nan’, ‘1.#IND’, ‘1.#QNAN’, ‘<NA>’, 
  ‘N/A’, ‘NA’, ‘NULL’, ‘NaN’, ‘n/a’, ‘nan’, ‘null’.


Internal datatypes
------------------

Acknowledge and cite Annotable
------------------------------

.. note:: We follow the heuristic that peptides are on the rows, 
    and samples are on the columns.
    This is not for any great reason other than we've seen this done 
    `historically <TODO>`_,
    and often we're performing vectorized operations on the samples, making it slightly faster 
    to make the enrichment matrix more "tall" than "wide". This may change in the future.

