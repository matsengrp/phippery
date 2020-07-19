"""
@File: analysis.py

@Author: Jared Galloway
"""


# dependencies
import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

# built-in
import itertools


def make_data_loader_infinite(data_loader):
    """With this we can always just ask for more data with next(), going
    through minibatches as guided by DataLoader."""
    for loader in itertools.repeat(data_loader):
        for data in loader:
            yield data


class Analysis:
    """A wrapper class for training models."""

    def __init__(
        self,
        model,
        train_data_list,
        val_data,
        batch_size=1,
        learning_rate=5e-3,
        device="cpu",
    ):
        self.batch_size = batch_size
        self.device = torch.device(device)
        self.learning_rate = learning_rate
        self.model = model
        self.val_data = val_data
        self.train_datasets = train_data_list
        self.train_loaders = [
            DataLoader(train_dataset, batch_size=self.batch_size, shuffle=True)
            for train_dataset in train_data_list
        ]
        self.train_infinite_loaders = [
            make_data_loader_infinite(train_loader)
            for train_loader in self.train_loaders
        ]

    def train(
        self, epoch_count, loss_fn, patience=10, min_lr=1e-5, loss_weight_span=None
    ):
        batch_loss, batch_val_loss = [], []
        """Train self.model using all the bells and whistles."""
        assert len(self.train_datasets) > 0
        batch_count = 1 + max(map(len, self.train_datasets)) // self.batch_size
        self.model.train()  # Sets model to training mode.
        self.model.to(self.device)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)

        def step_model():
            for i in range(batch_count):
                optimizer.zero_grad()
                per_batch_loss = 0.0
                for train_infinite_loader in self.train_infinite_loaders:

                    batch = next(train_infinite_loader)
                    samples = batch["samples"].to(self.device)
                    predictions = self.model(samples.float())
                    targets = batch["targets"].to(self.device)
                    loss = loss_fn(predictions, targets)
                    per_batch_loss += loss.item()

                    # Note that here we are using gradient accumulation: calling
                    # backward for each loader before clearing the gradient via
                    # zero_grad. See, e.g. https://link.medium.com/wem03OhPH5
                    loss.backward()
                optimizer.step()

            val_samples = torch.tensor(self.val_data.samples).to(self.device)
            val_targets = torch.tensor(self.val_data.targets).to(self.device)
            val_predictions = self.model(val_samples.float())
            val_loss = loss_fn(val_predictions, val_targets).item()
            # TODO impliment early stop based on validation loss.

            return per_batch_loss, val_loss

        for e in range(epoch_count):
            loss, val_loss = step_model()
            batch_loss.append(loss)
            batch_val_loss.append(val_loss)
            print(f"on epoch: {e} of {epoch_count}. Loss {loss} - val loss: {val_loss}")

        return batch_loss, batch_val_loss
