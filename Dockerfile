
FROM ubuntu:22.04

ENV LANG=en_US.UTF-8

COPY setup.py     /opt/primer_finder/setup.py
COPY cfeproviral /opt/primer_finder/cfeproviral/
COPY scripts      /opt/primer_finder/scripts/

RUN sh -- /opt/primer_finder/scripts/install.sh

WORKDIR /w

ENTRYPOINT ["cfeproviral"]
