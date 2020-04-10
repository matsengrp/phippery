# TODO
class PhipData(object):
    """
    A class to organize, query, and transform
    PhIP-seq data.
    """

    def __init__(
        self,
        counts_tables,  # 3D array-like
        peptide_metadata,  # 2D array-like
        sample_metadata,
    ):  # 2D array-like

        self.counts_tables = counts_tables
        self.peptide_metadata = peptide_metadata
        self.sample_metadata = sample_metadata
        # TODO
        # 1. Do some basic input format checking
        # self._assert_data_consistancy(counts_tables, peptide_metadata, sample_metadata)

        # 2. Infer Mock IP's and library control

        # 3. ????

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
