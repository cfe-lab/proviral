import os
import sys
import pandas as pd
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

import utils
from file import File
from fasta import Fasta


def test_get_softclip_start():
    example_aln = cwd / 'softclip.sam'
    example_aln = utils.load_samfile(example_aln)
    size, op = example_aln.iloc[0]['cigar'][0]
    size = int(size)
    query = Fasta(cwd / 'query.fasta')
    sequence = None
    for h, s in query:
        sequence = s
    assert sequence[:size] == utils.mod_hxb2[:size]
    assert sequence[size:] == utils.mod_hxb2[-len(sequence) + size:]