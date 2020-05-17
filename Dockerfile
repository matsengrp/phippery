FROM continuumio/anaconda:2019.10

RUN apt-get update -q && \
    apt-get install -y -q --no-install-recommends \
        build-essential

COPY environment.yml .

RUN /opt/conda/bin/conda env create -f environment.yml
