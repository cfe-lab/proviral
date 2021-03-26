import os
import pytest
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

from gene_splicer.probe_finder import unpack_mixtures_and_reverse


@pytest.fixture
def setup():
    return {'cwd': Path(os.path.realpath(__file__)).parent}


# Personally I had trouble with the unpack_mixtures_and_reverse() taking an inordinate amount of time to complete
def test_intractable(setup):
    seq = None
    with open(setup['cwd'] / 'data' / 'intractable_sequence.txt') as f:
        seq = f.read().strip()
    unpack_mixtures_and_reverse(seq)
