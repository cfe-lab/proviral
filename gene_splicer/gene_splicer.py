import utils
import argparse
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description='Splice out HIV genes from a sequence',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query_fasta', help='A query sequence')
    # parser.add_argument('target_fasta', help='The reference sequence')
    parser.add_argument('samfile', help='The sam alignment of query to target')
    # parser.add_argument('annotation', help='A CSV file with columns: "gene, start, stop" defining gene regions with respect to the target')
    # parser.add_argument('-o',
    #                     '--outpath',
    #                     help='The path to save the output',
    #                     default=os.getcwd())
    return parser.parse_args()


def main():
    args = parse_args()
    query = next(utils.read_fasta(args.query_fasta))[1]
    target = utils.mod_hxb2
    samfile = utils.load_samfile(args.samfile)
    results = utils.splice_genes(query, target, samfile, utils.mod_annot)
    genes = utils.coords_to_genes(results, query)
    genes_path = Path(args.query_fasta).resolve().parent / 'genes.fasta'
    with open(genes_path, 'w') as o:
        for gene, seq in genes.items():
            o.write(f'>{gene}\n{seq}\n')


if __name__ == '__main__':
    main()