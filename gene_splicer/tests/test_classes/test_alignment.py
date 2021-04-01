import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

import gene_splicer.utils as utils
from fasta import Fasta
from gene_splicer.alignment import Alignment


def test_alignment(tmp_path):
    data = cwd.parent / 'data'
    example = data / 'example2'
    utils.load_samfile(example / 'alignment.sam')
    # Read in the target sequence
    target = None
    for header, seq in Fasta(example / 'target.fasta'):
        target = seq
    # Read in the query sequence
    query = None
    for header, seq in Fasta(example / 'query.fasta'):
        query = seq
    aln = Alignment(target, query, tmp_path)
    test_aln = utils.load_samfile(aln.path)
    valid_aln = utils.load_samfile(example / 'alignment.sam')
    assert test_aln.equals(valid_aln)
