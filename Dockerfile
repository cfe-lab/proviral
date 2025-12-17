
FROM ubuntu:22.04

ENV LANG=en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

COPY pyproject.toml /opt/cfeproviral/pyproject.toml
COPY cfeproviral /opt/cfeproviral/cfeproviral/
COPY scripts      /opt/cfeproviral/scripts/

RUN sh -- /opt/cfeproviral/scripts/install.sh

WORKDIR /w

ENTRYPOINT ["cfeproviral"]
