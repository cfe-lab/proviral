#! /bin/sh

set -e

echo ===== Installing hivseqinr ===== >/dev/null

set -x

uvdo run -- cfeproviral hivseqinr /opt/hivseqinr

