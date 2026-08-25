#!/bin/sh
# Build a Singularity image from the freshly built docker image, the same way
# MiCall does: docker save -> definition file with Bootstrap: docker-archive ->
# singularity build. Results land in ./simgs with a proviral-latest.sif symlink.
set -eu

cd -- "$(git rev-parse --show-toplevel)"
mkdir -p simgs

image_name="$(sh scripts/docker_build.sh)"
container_sha="$(docker inspect --format '{{.Id}}' "$image_name" | sed 's/^sha256://; s/^\(.\{12\}\).*/\1/')"
archive_path="simgs/proviral-$container_sha.tar"
definition_path="simgs/proviral-$container_sha.def"
image_path="simgs/proviral-$container_sha.sif"

echo "Saving docker archive to $archive_path"
docker save --output "$archive_path" "$image_name"

cat > "$definition_path" <<EOF
Bootstrap: docker-archive
From: ./$archive_path

%help
    Search proviral consensus sequences for primers, then use CFEIntact to
    decide if the genomes are complete.

    This Singularity container can be run on Kive: http://cfe-lab.github.io/Kive


%labels
    MAINTAINER BC CfE in HIV/AIDS https://github.com/cfe-lab/
    KIVE_INPUTS sample_info_csv contigs_csv conseqs_csv cascade_csv
    KIVE_OUTPUTS outcome_summary_csv conseqs_primers_csv contigs_primers_csv \\
        table_precursor_csv proviral_landscape_csv detailed_results_tar
    KIVE_THREADS 1
    KIVE_MEMORY 6000

%environment
    export LANG=en_US.UTF-8

%runscript
    cd -- /w
    cfeproviral sample --cfeintact "\$@"
EOF

echo "Building $image_path"
singularity build "$image_path" "$definition_path"
ln -sf "$(basename "$image_path")" simgs/proviral-latest.sif
echo "Built $image_path"
