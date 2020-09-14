import os
import sys
from pathlib import Path

from pandas.core.arrays.period import validate_dtype_freq

cwd = Path(os.path.realpath(__file__)).parent
data = cwd.parent / 'data'
sys.path.append(str(cwd.parent.parent))

import utils
from fasta import Fasta
from alignment import Alignment


def test_alignment():
    example = data / 'example2'
    valid_aln = utils.load_samfile(example / 'alignment.sam')
    # Read in the target sequence
    target = None
    for header, seq in Fasta(example / 'target.fasta'):
        target = seq
    # Read in the query sequence
    query = None
    for header, seq in Fasta(example / 'query.fasta'):
        query = seq
    aln = Alignment(target, query, cwd / 'tmp')
    test_aln = utils.load_samfile(aln.path)
    valid_aln = utils.load_samfile(example / 'alignment.sam')
    assert test_aln == valid_aln