
FROM ubuntu:22.04

ENV LANG=en_US.UTF-8

COPY setup.py     /opt/primer_finder/setup.py
COPY gene_splicer /opt/primer_finder/gene_splicer/
COPY scripts      /opt/primer_finder/scripts/

RUN sh -- /opt/primer_finder/scripts/install.sh

ENTRYPOINT ["gene_splicer_sample", "--hivseqinr", "/opt/hivseqinr"]
