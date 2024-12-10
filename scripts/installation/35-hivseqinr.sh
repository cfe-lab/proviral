#! /bin/sh

set -e

echo ===== Installing hivseqinr ===== >/dev/null

set -x

export DEBIAN_FRONTEND=noninteractive

apt-get install -y libz-dev libcurl4-openssl-dev libxml2-dev
apt-get install --no-install-recommends -y r-base
Rscript /opt/primer_finder/gene_splicer/configure_r.sh
python3 -m gene_splicer.hivseqinr /opt/hivseqinr
