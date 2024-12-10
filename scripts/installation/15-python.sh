#! /bin/sh

set -e

echo ===== Installing Python ===== >/dev/null

set -x

apt-get install -y python3 python3-pip
