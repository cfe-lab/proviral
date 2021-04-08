import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

from fasta import Fasta
import gene_splicer.utils as utils


def test_fasta():
    fasta = Fasta()
    expected_header = '>mod_hxb2'
    expected_seq = 'GGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGC'
    fasta.write(expected_header + '\n')
    fasta.write(expected_seq + '\n')
    fasta.seek(0)
    lines = fasta.readlines()
    assert lines == [expected_header + '\n', expected_seq + '\n']
    fasta.seek(0)
    for header, seq in fasta:
        assert header == expected_header
        assert seq == expected_seq
