import os
import io
import sys
import csv
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

import gene_splicer.utils as utils


def test_merge_coords():
    coords1 = {'A': [10, 20], 'B': [30, 40]}
    coords2 = {'A': [10, 30], 'B': [30, 40]}
    merged1 = utils.merge_coords(coords2, coords1)
    merged2 = utils.merge_coords(coords1, coords2)
    assert merged1 == merged2
    assert merged1 == coords2


def test_getSamplesFromCascade():
    # Write an in-memory cascade file
    cascade = io.StringIO(newline='')
    fieldnames = (
        'sample',
        'demultiplexed',
        'v3loop',
        'g2p',
        'prelim_map',
        'remap',
        'aligned',
    )
    writer = csv.DictWriter(cascade, fieldnames)
    writer.writeheader()
    for i in range(10):
        writer.writerow({
            'sample': i,
            'demultiplexed': i,
            'v3loop': i,
            'g2p': i,
            'prelim_map': i,
            'remap': i,
            'aligned': i
        })
    cascade.seek(0)

    # Test the function
    samples = utils.getSamplesFromCascade(cascade)
    assert len(samples) == 10
    for i in range(10):
        assert samples[str(i)] == i