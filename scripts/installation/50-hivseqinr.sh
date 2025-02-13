#! /bin/sh

set -e

echo ===== Installing hivseqinr ===== >/dev/null

set -x

python3 -m cfeproviral.hivseqinr /opt/hivseqinr

