#! /bin/sh

set -e

echo ===== Installing R dependencies needed for hivseqinr ===== >/dev/null

set -x

export DEBIAN_FRONTEND=noninteractive

apt-get install -y libz-dev libcurl4-openssl-dev libxml2-dev
apt-get install --no-install-recommends -y r-base
Rscript /opt/primer_finder/cfeproviral/configure_r.sh
