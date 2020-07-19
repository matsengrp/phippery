"""
@File: loss.py

@Author: Jared Galloway
"""

from torch import nn


def mse(y_true, y_predicted):
    return nn.functional.mse_loss(y_true.squeeze(), y_predicted.squeeze())
