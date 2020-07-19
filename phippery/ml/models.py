"""
@File: models.py

@Author: Jared Galloway

A file to collect some machine learning model
architectures.
"""


# dependencies
from torch import nn


class Net(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(Net, self).__init__()  # Inherited from the parent class nn.Module
        self.fc1 = nn.Linear(
            input_size, hidden_size
        )  # 1st Full-Connected Layer: 784 (input data) -> 500 (hidden node)
        self.relu = nn.ReLU()  # Non-Linear ReLU Layer: max(0,x)
        self.fc2 = nn.Linear(
            hidden_size, num_classes
        )  # 2nd Full-Connected Layer: 500 (hidden node) -> 10 (output class)

    def forward(self, x):  # Forward pass: stacking each layer together
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out


class DropoutNet(nn.Module):
    def __init__(self, input_size, hidden_size, num_classes):
        super(DropoutNet, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.drop1 = nn.Dropout(0.1)
        self.s1 = nn.Sigmoid()
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.drop2 = nn.Dropout(0.1)
        self.s2 = nn.Sigmoid()
        self.fc3 = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        out = self.fc1(x)
        out = self.drop1(out)
        out = self.s1(out)
        out = self.fc2(out)
        out = self.drop2(out)
        out = self.s2(out)
        out = self.fc3(out)
        return out
