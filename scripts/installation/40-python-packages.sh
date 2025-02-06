#! /bin/sh

set -e

echo ===== Installing Python packages ===== >/dev/null

set -x

pip3 install /opt/primer_finder
python3 -m gene_splicer.hivseqinr /opt/hivseqinr
