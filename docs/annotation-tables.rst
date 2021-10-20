
.. _sec_anno_intro:

=================
Annotation Tables
=================

.. _sec_pipeline_anno:

Informative annotations for the enrichment matrix rows (peptides), 
and columns (samples) are at the core of our approach to analysis of this
type of data. The high-throughput, brute-force nature of this protocol
results in what may be most clearly understand as *many* small experiments 
- each giving us detailed information about binding affinity 
of sampled antibodies across a vast combinations of viral (or other) proteomes.

Because this information is so useful for computing various transformations on the data,
we require the user provides these annotations tables in a strait-forward, 
albeit *specific* format. It should be noted that there are really no r
equired annotations to use a small subset of `phippery's` useful functions, 
but the more information you provide outlining your specific dataset, the 
powerful the software becomes. Here, we outline the formatting requirements, 
as well as some useful annotations you might consider having prepared before 
using the tools presented here. 

.. note:: While `phippery` has no required annotation features for either peptides, 
    or samples, the `Nextflow pipeline <TODO>`_, *does* require a 
    few special columns necessary for performing the alignment steps correctly. 
    More on this, below.

++++++++++++++++++++++++
Annotation Table Gotchas
++++++++++++++++++++++++

.. note:: We follow the heuristic that peptides are on the rows, 
    and samples are on the columns.
    This is not for any great reason other than we've seen this done 
    `historically <TODO>`_,
    and often we're performing vectorized operations on the samples, making it slightly faster 
    to make the enrichment matrix more "tall" than "wide". This may change in the future.

.. note:: When dealing with missing values in the annotation tables, we use the 
    `pd.convert_dtypes <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.convert_dtypes.html>`_
    function to best allow for missing annotations, while maintaining the integrety of
    the inferred datatype. It is highly reccomended you stay consistant with datatypes for feature annotations,
    i.e. try not to mix values like `1` (integer), `6.7` (float), and `hello_world` (string) in any one of the columns.
  
    for missing data of any type, The following values will be interpretted as `NaN`; ‘’, ‘#N/A’, ‘#N/A N/A’, 
    ‘#NA’, ‘-1.#IND’, ‘-1.#QNAN’, ‘-NaN’, ‘-nan’, ‘1.#IND’, ‘1.#QNAN’, ‘<NA>’, 
    ‘N/A’, ‘NA’, ‘NULL’, ‘NaN’, ‘n/a’, ‘nan’, ‘null’.

