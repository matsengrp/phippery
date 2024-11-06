FROM quay.io/hdc-workflows/ubuntu:20.04

# bust cache
ADD http://date.jsontest.com /etc/builddate

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
      version "1.3.0" \
      description "Common PhIP-Seq Workflows"

# install needed tools
RUN apt-get update --fix-missing -qq && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y -q \
    git \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    python3 \
    python3-pip \ 
    python3-venv \
    zip \
    wget

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# install phippery
RUN pip install git+https://github.com/matsengrp/phippery@1.3.0

# install pre-build binary Bowtie1.3
RUN curl -fksSL https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.1-linux-x86_64.zip \
    --output bowtie-1.3.0-linux-x86_64.zip \
    && unzip bowtie-1.3.0-linux-x86_64.zip \
    && (cd /usr/bin/ && ln -s /bowtie-1.3.0-linux-x86_64/* ./)


# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.0/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.0 && \
    make all all-htslib && make install install-htslib

