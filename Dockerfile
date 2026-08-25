FROM ubuntu:24.04

ENV LANG=en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

# Each installation step runs in its own layer, so rebuilt images reuse the
# cached layers of every step that didn't change.
COPY scripts/installation/ /opt/cfeproviral/scripts/installation/

RUN sh -- /opt/cfeproviral/scripts/installation/10-prerequisites.sh

RUN sh -- /opt/cfeproviral/scripts/installation/20-blast.sh
RUN sh -- /opt/cfeproviral/scripts/installation/20-mafft.sh
RUN sh -- /opt/cfeproviral/scripts/installation/30-minimap2.sh

COPY cfeproviral/configure_r.sh /opt/cfeproviral/cfeproviral/configure_r.sh
RUN sh -- /opt/cfeproviral/scripts/installation/35-r-packages.sh

COPY pyproject.toml cfeproviral/ /opt/cfeproviral/
RUN sh -- /opt/cfeproviral/scripts/installation/40-python-packages.sh

RUN sh -- /opt/cfeproviral/scripts/installation/50-hivseqinr.sh

RUN sh -- /opt/cfeproviral/scripts/installation/90-cleanup.sh

WORKDIR /w

ENTRYPOINT ["cfeproviral"]
