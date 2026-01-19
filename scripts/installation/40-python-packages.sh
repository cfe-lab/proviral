#! /bin/sh

set -e

echo ===== Installing Python packages ===== >/dev/null

set -x

pip3 install --break-system-packages /opt/cfeproviral
cfeproviral --version
