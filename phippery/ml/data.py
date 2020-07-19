"""
@File: data.py

@Author: Jared Galloway
"""


# dependencies
import xarray as xr
from torch.utils.data import Dataset

# built-in
import random
import numpy as np

# local
from phippery.utils import get_all_sample_metadata_factors
from phippery.utils import iter_sample_groups


class PhipDataset(Dataset):
    """
    A very simple dataset class
    """

    def __init__(self, samples, targets):
        self.samples = samples
        self.targets = targets

    def __getitem__(self, idxs):
        return {"samples": self.samples[idxs], "targets": self.targets[idxs]}

    def __len__(self):
        return self.samples.shape[0]


def get_train_test_split_by(
    ds, sample_target="patient_status", percent_test=0.1, bio_id_column="sample_ID"
):
    """
    choose a target from sample metadata and create a train/test split
    based on enrichments from the dataset.
    """

    if sample_target not in ds.sample_metadata.values:
        raise ValueError(f"{sample_target} is not in sample_metadata coordinate")

    train_pds_dict = {}
    test_X, test_Y = [], []
    for label, (patient_status, group) in enumerate(
        iter_sample_groups(ds, sample_target)
    ):
        unique_bio = get_all_sample_metadata_factors(group, bio_id_column)
        n_test = int(len(unique_bio) * percent_test)
        test_bio = random.sample(unique_bio, n_test)
        train_X, train_Y = [], []
        for bio_replicate, bio_group in iter_sample_groups(group, bio_id_column):
            if bio_replicate in test_bio:
                for replicate in bio_group.sample_id.values:
                    test_X.append(bio_group.counts.loc[:, replicate])
                    test_Y.append(label)
            else:
                for replicate in bio_group.sample_id.values:
                    train_X.append(bio_group.counts.loc[:, replicate])
                    train_Y.append(label)

        train_pds_dict[patient_status] = PhipDataset(
            np.array(train_X), np.array(train_Y)
        )
    test_pds = PhipDataset(np.array(test_X), np.array(test_Y))
    return train_pds_dict, test_pds
