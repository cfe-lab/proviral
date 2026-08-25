FROM ubuntu:24.04

ENV LANG=en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Prerequisites.
RUN mkdir -p /w \
    && apt-get update -qq \
    && apt-get install -y build-essential unzip git wget \
        fontconfig libbz2-dev liblzma-dev libssl-dev \
        libffi-dev libsqlite3-dev tar python3 python3-pip

RUN apt-get install -y ncbi-blast+

# MAFFT is an aligner program, like clustalw, mappy. It is used by CFEIntact.
RUN apt-get install -y mafft

RUN apt-get install -y minimap2

# R dependencies needed for hivseqinr.
COPY cfeproviral/configure_r.sh /opt/cfeproviral/cfeproviral/configure_r.sh
RUN apt-get install -y libz-dev libcurl4-openssl-dev libxml2-dev \
    && apt-get install --no-install-recommends -y r-base \
    && Rscript /opt/cfeproviral/cfeproviral/configure_r.sh

COPY pyproject.toml cfeproviral/ /opt/cfeproviral/
RUN pip3 install --break-system-packages /opt/cfeproviral \
    && cfeproviral --version

RUN cfeproviral hivseqinr /opt/hivseqinr

RUN apt-get remove -y wget git build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /w

ENTRYPOINT ["cfeproviral"]
