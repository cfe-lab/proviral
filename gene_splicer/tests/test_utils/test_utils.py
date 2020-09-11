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
    query_fasta = Fasta(cwd / 'query.fasta')
    query_sequence = None
    for h, s in query_fasta:
        query_sequence = s
    assert query_sequence[:size] == utils.mod_hxb2[:size]
    assert query_sequence[size:] == utils.mod_hxb2[-len(query_sequence) +
                                                   size:]
    # softclip_start = utils.get_softclip_start(query_sequence, example_aln)
    # expected_softclip_start = query_sequence[:size]
    # assert softclip_start == expected_softclip_start
    outpath = cwd / 'deleteme'
    # softclip_start = utils.get_softclip_start(utils.mod_hxb2, query_sequence,
    #                                           example_aln, outpath)
    softclip_start = utils.align(utils.mod_hxb2, query_sequence, outpath)
