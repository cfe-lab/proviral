import os
import pytest
from pathlib import Path

from cfeproviral.probe_finder import ProbeFinder


@pytest.fixture
def currentdir():
    cwd = Path(os.path.realpath(__file__)).parent
    return cwd


def test_probe_finder(currentdir):
    # Try to find TCGA, with 10 "X"s on either side
    target = 'CCCCCCCCCCTCGACCCCCCCCCC'
    query = 'TCGA'
    finder = ProbeFinder(query, target)
    assert finder.contig_match == 'TCGA'
