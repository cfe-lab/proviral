import utils
import argparse
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Splice out HIV genes from a sequence',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query_fasta', help='A fasta file containing one or more query sequences')
    return parser.parse_args()


def run(query_fasta, outdir=Path(os.getcwd()).resolve()):
    target = utils.mod_hxb2
    for query_name, query_seq in utils.read_fasta(query_fasta):
        # Splitting by '::' is quite specific, make sure primer_finder joins using this
        samfile_path = utils.align(target, query_seq, query_name[1:].split('::')[1], outdir=outdir)
        samfile = utils.load_samfile(samfile_path)
        results = utils.splice_genes(query_seq, target, samfile, utils.mod_annot)
        genes = utils.coords_to_genes(results, query_seq)
        genes_path = samfile_path.parent / 'genes.fasta'
        with open(genes_path, 'w') as o:
            for gene, seq in genes.items():
                o.write(f'>{gene}\n{seq}\n')


def main():
    args = parse_args()
    run(args.query_fasta)


if __name__ == '__main__':
    main()