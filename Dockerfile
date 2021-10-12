FROM python:3.7
RUN python -m pip install --upgrade pip
COPY . .
RUN python setup.py install
