#! /bin/sh

set -e

echo ===== "Installing python's uv" ===== >/dev/null

set -x

wget -q https://astral.sh/uv/install.sh -O "/tmp/uv-install.sh"
sh "/tmp/uv-install.sh"

cp -v -- ~/.local/bin/uv ~/.local/bin/uvx /bin
uv --version
