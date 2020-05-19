FROM continuumio/anaconda:2019.10

RUN apt-get update -q && \
    apt-get install -y -q --no-install-recommends \
        build-essential

COPY . .
RUN conda install --quiet --yes --file requirements2.txt --channel conda-forge && conda clean --all -f -y
RUN pip install -e .
