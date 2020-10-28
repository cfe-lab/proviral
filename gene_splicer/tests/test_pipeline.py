import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent))

import primer_finder


def test_v3_removal():
    conseqs = open(cwd / 'data' / 'v3only' / 'conseq.csv')
    outpath = cwd
    run_name = 'testrunv3'
    outfile = primer_finder.find_primers(conseqs, outpath, run_name)
    with open(outfile) as output:
        lines = output.readlines()
    assert len(lines) == 2
    assert lines[1].split(',')[3] == 'is V3 sequence'