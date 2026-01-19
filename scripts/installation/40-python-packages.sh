#! /bin/sh

set -e

echo ===== Installing Python packages ===== >/dev/null

set -x

uvdo sync
uvdo run -- cfeproviral --version
