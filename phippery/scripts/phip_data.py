"""
file: phip_data.py

@Author: Jared Galloway
some inspiration: https://github.com/lasersonlab/phip-stat

This file contains the code
to take in peptide enrichent data in the form of
counts files, merges them, and then checks consistancy

(first column is the peptide id and
second is the alignment counts)
"""


# TODO
class PhipData(object):
    """
    A class to organize, query, and transform
    PhIP-seq data.


    """

    def __init__(counts_dir):

        pass

    def _assert_data_consistancy_stub(self, counts, sample_md, peptide_md):
        """
        'private' method to do some sanity checking and error handling
        on the inputs provided. This will provide some
        useful warning error messages, as well as provide a platform for
        testing that the input will be valid for all methods below.
        """
        pass

    def get_mock_ip_stub(self):
        """
        """
        pass

    def get_library_control_stub(self):
        """
        """
        pass

    def dump_stub(self, output_path):
        """
        (possibly) Organize the three tables into a dictionary-like structure
        Use pickle/numpy/pandas to dump all the data
        and make it easy to unpack and load again into the same object.
        """
        pass

    def technical_replicate_correlations_stub(self):
        """
        compute and return all correlations between
        available technical replicates.
        """
        pass
