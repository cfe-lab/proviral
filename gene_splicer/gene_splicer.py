import argparse
import os
from pathlib import Path
import gene_splicer.utils as utils


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
        alignment_path = utils.align(target,
                                     query_seq,
                                     query_name[1:].split('::')[1],
                                     outdir=outdir)
        if not alignment_path:
            print(f'Could not align {query_name}, aligner not available')
            continue
        samfile = utils.load_samfile(alignment_path)
        coords = utils.splice_genes(query_seq, target, samfile,
                                    utils.mod_annot)
        # Try to get softclipped region if there is one
        softclipped_coords = utils.sequence_to_coords(query_seq, target,
                                                      alignment_path,
                                                      utils.mod_annot)
        if softclipped_coords:
            coords = utils.merge_coords(softclipped_coords, coords)
        genes = utils.coords_to_genes(coords, query_seq)
        genes_path = alignment_path.parent / 'genes.fasta'
        with open(genes_path, 'w') as o:
            for gene, seq in genes.items():
                o.write(f'>{gene}\n{seq}\n')


def main():
    args = parse_args()
    run(query_fasta=args.query_fasta, outdir=args.outdir)


if __name__ == '__main__':
    main()
