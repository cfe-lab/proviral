#!/bin/sh
# Build the proviral docker image and tag it for GHCR using the version
# reported by git describe (leading 'v' stripped), e.g.
# ghcr.io/cfe-lab/proviral:2.5.3-81-gb07fdfc
set -eu

cd -- "$(git rev-parse --show-toplevel)"
version="$(git describe | sed 's/^v//')"
image_name="ghcr.io/cfe-lab/proviral:$version"
docker build --tag "$image_name" .
echo "$image_name"
