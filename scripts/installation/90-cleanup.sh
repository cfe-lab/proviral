#! /bin/sh

set -e

echo ===== Clean up ===== >/dev/null

set -x

apt-get remove -y wget build-essential
apt-get clean
rm -rf /var/lib/apt/lists/*
