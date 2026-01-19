
FROM ubuntu:22.04

ENV LANG=en_US.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

COPY . /opt/cfeproviral

RUN sh -- /opt/cfeproviral/scripts/install.sh

WORKDIR /w

ENTRYPOINT ["uvrun", "cfeproviral"]
