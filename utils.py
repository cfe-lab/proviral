import os
import re
from pathlib import Path
from csv import DictWriter, DictReader


def reverse_and_complement(seq):
    return ''.join(complement_dict[nuc] for nuc in reversed(seq))


def write_annot(adict, filename):
    """
        adict: {gene (str): [start (int), stop (int)], ...}
    """
    with open(filename, 'w', newline='') as o:
        fieldnames = ['gene', 'start', 'stop']
        writer = DictWriter(o, fieldnames)
        writer.writeheader()
        for gene, (start, stop) in adict.items():
            writer.writerow({
                'gene': gene.upper(),
                'start': start,
                'stop': stop
            })


def write_fasta(adict, filename):
    with open(filename, 'w') as o:
        for key, value in adict.items():
            o.write(f'>{key}\n{value}\n')


def read_csv(csvfile):
    with open(csvfile) as f:
        reader = DictReader(f)
        for row in reader:
            yield row


def csv_to_bed(csvfile, target_name='HXB2', offset_start=0, offset_stop=0):
    with open(csvfile) as f:
        reader = DictReader(f)
        with open(f'{csvfile}.bed', 'w', newline='') as o:
            fieldnames = [
                'chrom',
                'chromStart',
                'chromEnd',
                'name',
                'score',
                'strand',
                'thickStart',
                'thickEnd'
            ]
            writer = DictWriter(o, fieldnames, delimiter='\t')
            for row in reader:
                writer.writerow({
                    'chrom': target_name,
                    'chromStart': int(row['start']) + offset_start,
                    'chromEnd': int(row['stop']) + offset_stop,
                    'name': row['gene'].upper(),
                    'score': 0,
                    'strand': '+',
                    'thickStart': int(row['start']) + offset_start,
                    'thickEnd': int(row['stop']) + offset_stop
                })


def split_cigar(row):
    pattern = re.compile('(\d+)([A-Z])')
    cigar = re.findall(pattern, row[5])
    return cigar


def read_fasta(fasta_file):
    with open(fasta_file) as fasta:
        name, seq = None, []
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            yield (name, ''.join(seq))


def get_genes(annot, pos, offset_start=0, offset_end=0, offset=None):
    genes = []
    if offset is not None:
        offset_start = offset_end = offset

    """
        annot (dict) => {gene: [start, end], ...}
        ref (str),
        pos (int)
    """
    for gene in annot:
        if annot[gene][0] - offset_start <= pos <= annot[gene][1] - offset_end:
            genes.append(gene)
    return genes

# Define some global variables
mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}

complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'S', 'R': 'Y', 'K': 'M', 'Y': 'R', 'S': 'W', 'M': 'K',
                   'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
                   '*': '*', 'N': 'N', '-': '-'}

cwd = Path(os.path.realpath(__file__)).parent

hxb2 = next(read_fasta(cwd / 'hxb2.fasta'))[1]