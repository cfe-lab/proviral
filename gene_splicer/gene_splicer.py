import utils
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description='Splice out HIV genes from a sequence',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('query_fasta', help='A query sequence')
    parser.add_argument('target_fasta', help='The reference sequence')
    parser.add_argument('annotation', help='A CSV file with columns: "gene, start, stop" defining gene regions with respect to the target')
    parser.add_argument('-o',
                        '--outpath',
                        help='The path to save the output',
                        default=os.getcwd())
    return parser.parse_args()
