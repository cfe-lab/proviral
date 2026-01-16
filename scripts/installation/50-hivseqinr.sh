#! /bin/sh

set -e

echo ===== Installing hivseqinr ===== >/dev/null

set -x

uv --project /opt/cfeproviral run -- cfeproviral hivseqinr /opt/hivseqinr

