import pytest
import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

import gene_splicer.utils as utils
from gene_splicer.probe_finder import ProbeFinder
from gene_splicer.probe_finder import TARGET_SEQUENCES
from gene_splicer.primer_finder import primers, validate_primer


@pytest.fixture
def data():
    cwd = Path(os.path.realpath(__file__)).parent
    seq = utils.load_yaml(cwd / 'data' / '3prime_end.yaml')
    return {'cwd': cwd, 'seq': seq}


def test_3prime_end(data):
    # finder = ProbeFinder(data['seq'], TARGET_SEQUENCES['round2_rev_primer'])
    hxb2_target_start = primers['rev']['hxb2_start'] - 100
    hxb2_target_end = primers['rev']['hxb2_end']
    hxb2_target_seq = utils.hxb2[hxb2_target_start:hxb2_target_end]
    finder = ProbeFinder(hxb2_target_seq, data['seq'])
    finder.start += hxb2_target_start
    print()
    print(finder)
    validated = validate_primer(finder, data['seq'], primers['rev'])
    print()
    print(validated)