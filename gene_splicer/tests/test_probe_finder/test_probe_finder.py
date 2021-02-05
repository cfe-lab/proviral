import os
import pytest
from pathlib import Path

from gene_splicer.probe_finder import ProbeFinder
import gene_splicer.utils as utils
from gene_splicer.fasta import Fasta


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


def test_real_example1(currentdir):
    example = currentdir / 'ignore_example1'
    fasta = Fasta(example / 'query.fasta')
    for h, s in fasta:
        query = s
    target = utils.mod_hxb2
    aln_path = example / 'alignment.sam'
    coords = utils.sequence_to_coords(query, target, aln_path, utils.mod_annot)
    genes = utils.coords_to_genes(coords, query)
    assert 'x1' in genes
    assert genes['x1'] == 'AGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG'
