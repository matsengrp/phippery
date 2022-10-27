FROM quay.io/hdc-workflows/ubuntu:20.04

# bust cache
ADD http://date.jsontest.com /etc/builddate

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
    version "1.1.2" \
    description "Common PhIP-Seq Workflows"

# install ubuntu basics
RUN apt-get update --fix-missing -qq && \
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
ADD . /phippery
RUN (cd /phippery && pip install -e ".[dev]")
