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


def main():
    args = parse_args()
    target = utils.mod_hxb2
    for query_name, query_seq in utils.read_fasta(args.query_fasta):
        samfile = utils.align(target, query_seq, query_name[1:])
        results = utils.splice_genes(query_seq, target, samfile, utils.mod_annot)
        genes = utils.coords_to_genes(results, query_seq)
        genes_path = samfile.parent / 'genes.fasta'
        with open(genes_path, 'w') as o:
            for gene, seq in genes.items():
                o.write(f'>{gene}\n{seq}\n')


if __name__ == '__main__':
    main()