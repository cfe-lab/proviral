#! /bin/sh

set -e

echo ===== Installing Prerequisites ===== >/dev/null

set -x

apt-get update -qq
apt-get install -y build-essential unzip git wget \
        fontconfig libbz2-dev liblzma-dev libssl-dev \
        libffi-dev libsqlite3-dev tar
