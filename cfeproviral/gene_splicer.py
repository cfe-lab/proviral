import argparse
import os
import sys
from pathlib import Path
import cfeproviral.utils as utils


def parse_args():
    parser = argparse.ArgumentParser(
        description='Splice out HIV genes from a sequence',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'query_fasta',
        help='A fasta file containing one or more query sequences')
    parser.add_argument('--name',
                        help='A name for your analysis',
                        default='default_name')
    parser.add_argument('--outdir',
                        type=Path,
                        default=Path(os.getcwd()).resolve(),
                        help='Path to output files')
    return parser.parse_args()


def run(query_fasta, outdir):
    target = utils.mod_hxb2
    for query_name, query_seq in utils.read_fasta(query_fasta):
        # Splitting by '::' is quite specific, make sure primer_finder joins using this
        cigar_hits = utils.align(target,
                                 query_seq,
                                 query_name[1:].split('::')[1],
                                 outdir=outdir)
        if not cigar_hits:
            print(f'Could not align {query_name}, no alignment results')
            continue
        coords = utils.splice_genes(query_seq, target, cigar_hits,
                                    utils.mod_annot)
        # Try to get softclipped region if there is one
        softclipped_coords = utils.sequence_to_coords(query_seq, target,
                                                      cigar_hits,
                                                      utils.mod_annot)
        if softclipped_coords:
            coords = utils.merge_coords(softclipped_coords, coords)
        genes = utils.coords_to_genes(coords, query_seq)
        genes_path = outdir / query_name[1:].split('::')[1] / 'genes.fasta'
        with open(genes_path, 'w') as o:
            for gene, seq in genes.items():
                o.write(f'>{gene}\n{seq}\n')


def main():
    args = parse_args()
    run(query_fasta=args.query_fasta, outdir=args.outdir)


if __name__ == '__main__':
    main()
