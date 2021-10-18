FROM quay.io/matsengrp/python3.7
RUN python -m pip install --upgrade pip
COPY . .
RUN python setup.py install
