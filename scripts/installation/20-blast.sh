#! /bin/sh

set -e

echo ===== Installing blast ===== >/dev/null

set -x

apt-get install -y ncbi-blast+
