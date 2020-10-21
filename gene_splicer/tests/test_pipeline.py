import os
import sys
from pathlib import Path
from typing import runtime_checkable

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent))

import primer_finder


def test_v3_removal():
    conseqs = cwd / 'data' / 'v3only' / 'conseq.csv'
    outpath = cwd
    run_name = 'testrunv3'
    outfile = primer_finder.find_primers(conseqs, outpath, run_name)
    with open(outfile) as output:
        lines = output.readlines()
    assert len(lines) == 1