#! /bin/sh

set -e

echo ===== Installing Python packages ===== >/dev/null

set -x

uv --project /opt/cfeproviral sync
