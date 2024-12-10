#! /bin/sh

set -e

CURRENT_DIR="${0%/*}"

set -x

for SCRIPT in "$CURRENT_DIR/installation/"*
do
    sh -- "$SCRIPT"
done
