import os
import sys
import yaml
from pathlib import Path
import io

cwd = Path(os.path.realpath(__file__)).parent
inputs = cwd / 'inputs'
outputs = cwd / 'outputs'

import gene_splicer.utils as utils
import gene_splicer.primer_finder as primer_finder
import gene_splicer.gene_splicer as gene_splicer

utils.write_fasta({'MOD_HXB2': utils.mod_hxb2}, inputs / 'mod_hxb2.fasta')


class Args:
    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)


def test_intact():
    # For this test I will simply use the full HXB2 reference as my sequence, this should guarantee it's intact-ness
    alignment = utils.load_samfile(inputs / 'alignment.sam')
    coords = utils.splice_genes(utils.mod_hxb2, utils.mod_hxb2, alignment,
                                utils.mod_annot)
    genes = utils.coords_to_genes(coords, utils.mod_hxb2)
    genes_file = io.StringIO()
    utils.write_fasta(genes, genes_file)
    genes_file.seek(0)

    valid_genes = utils.load_yaml(inputs / 'valid_genes.yaml')
    for gene, seq in genes.items():
        assert seq == valid_genes[gene]
    # Mock hivseqinr results
    hivseqinr_resultsfile = io.StringIO()
    hivseqinr_resultsfile.write('SEQID,MyVerdict\n')
    hivseqinr_resultsfile.write('fakeref::contig,Intact\n')
    hivseqinr_resultsfile.seek(0)
    # Mock filtered (qc-passed) file
    filtered_csvfile = io.StringIO()
    filtered_csvfile.write(','.join(('reference', 'seqtype', 'sequence',
                                     'seqlen')) + '\n')
    filtered_csvfile.write(','.join(
        ('fakeref', 'contig', utils.mod_hxb2, str(len(utils.mod_hxb2)))))
    filtered_csvfile.seek(0)

    table_precursor_file = io.StringIO()
    utils.generate_table_precursor_2(hivseqinr_resultsfile, filtered_csvfile,
                                     genes_file, table_precursor_file)


def test_pipeline_sample1():
    conseq_path = cwd / 'inputs' / 'sample1' / 'conseq.csv'
    contigs_path = cwd / 'inputs' / 'sample1' / 'contigs.csv'
    outpath = cwd / 'outputs' / 'sample1'
    fasta_paths = primer_finder.run(conseqs_csv=open(conseq_path),
                                    contigs_csv=open(contigs_path),
                                    outpath=(outpath),
                                    extended_size=1000)
    for fasta in fasta_paths:
        args = Args(query_fasta=fasta, name='sample1', outpath=outpath)
        gene_splicer.run(fasta, args)