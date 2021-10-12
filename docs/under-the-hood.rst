

.. _sec_python_intro:

==============
Under the Hood
==============

The primary data strucure resulting from PhIP-Seq experiments is an *enrichment matrix*, 
X, with i rows and j columns. 
Commonly, row index represents a peptide that is displayed on a phage,
and each column represents a sample that was mixed with the entire phage library. 
After sequencing and demultiplexing each sample, we align the reads to the 
oligonucleotide reference library to observe a
count of aligned reads to each peptide.

Outside of the enrichment matrix, each *sample* in an experiment as well as each *peptide*
in the phage library used have number of important annotations required when
performing analysis tasks like model fitting, normalizing, and differential selection.
Additionally, the comparison across groups of virus proteins and 
sample types is crucial in many experiments. For large sample size experiments, 
it can be difficult to cross reference each of these groups before and
after analysis. 


Here, we take advantage of the powerful 
`xarray <http://xarray.pydata.org/en/stable/index.html>`_
approach to organizing all the Phip-Seq data along four primary coordinate 
dimensions which tie all sample/peptide enrichments to the respective annotations. 
Doing this allows us to store all the information without the error prone 
step of cross-checking separate dataframes, and without the
large storage scaling of using "Tall" dataframes.


