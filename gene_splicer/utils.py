import os
import re
from pathlib import Path
from csv import DictWriter, DictReader
from logger import logger


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


def modify_reference(refseq):
    newseq = refseq

    # Second round primers (regex for mixture positions)
    fwd = 'GCGCCCGAACAGGGAC.TGAAA.CGAAAG'
    rev = 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC'
    fwd_pos = re.search(fwd, newseq).start()
    rev_pos = [x for x in re.finditer(rev, newseq)][1].start()

    fwd_offset = fwd_pos + len(fwd)

    # Trim primers
    newseq = newseq[fwd_offset:rev_pos]

    assert newseq[:10] == 'GGAAACCAGA'

    # Reading frame correction (remove 1 nucleotide)
    newseq = newseq[:5107] + newseq[5108:]
    assert newseq[5102:5112] == 'CATTTCAGAA'

    # Fix premature stop codon
    pos = 8499
    assert newseq[pos] == 'A'
    newseq = newseq[:pos] + 'G' + newseq[pos+1:]
    assert newseq[pos] == 'G'

    return newseq


def modify_annot(annot):
    newannot = {}
    genes_of_interest = [
        'GAG',
        'POL',
        'PROTEASE',
        'RT',
        'INTEGRASE',
        'VIF',
        'APR',
        'ENV',
        'TAT_EXON1',
        'TAT_EXON2',
        'REV_EXON1',
        'REV_EXON2',
        'VPU',
        'VPR',
        'GP120',
        'GP41',
        'NEF',
        'PSI'
    ]
    for gene, (start, stop) in annot.items():
        if gene in genes_of_interest:
            newannot[gene] = [start, stop]
#     newannot = dict(annot)
    # Offset by second round fwd primer trim (666 trimmed from start)
    offset = -666
    for gene, (start, stop) in newannot.items():
        start += offset
        stop += offset
        if start <= 0:
            continue
        # Account for 1 nucleotide removal at index 5107 (pos 5108)
        if start >= 5108:
            start -= 1
            stop -= 1
        elif stop >= 5108:
            stop -=1
        # Convert to 0 base
        newannot[gene] = [start - 1, stop]

    return newannot


def splice_genes(query, target, samfile, annotation):
    results = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        for size, op in row['cigar']:
            size = int(size)
            print(f'size: {size}, op: {op}')
            print(f'target_pos: {target_pos}, query_pos: {query_pos}')
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                print('='*50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                print('='*50)
                continue
            elif op in ('M', '=', 'X'):
                for i in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        print(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    print(target_pos, query_pos, target_nuc, query_nuc, match, genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                        elif query_pos > results[gene][1]:
                                results[gene][1] = query_pos
                    query_pos += 1
                    target_pos += 1
            elif op == 'D':
                target_pos += size
                print('='*50)
                continue
            elif op == 'I':
                query_pos += size
                print('='*50)
                continue
            print('='*50)
        print('new alignment row'.center(50, '~'))
    return results


# This function has not been tested yet
def get_sequences(query, target, samfile, annotation):
    results = {}
    sequences = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        genes = get_genes(annotation, target_pos)
        query_pos = None
        for size, op in row['cigar']:
            print(f'size: {size}, op: {op}')
            print(f'target_pos: {target_pos}, query_pos: {query_pos}')
            size = int(size)
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                print('='*50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                print('='*50)
                continue
            elif op in ('M', '=', 'X'):
                for i in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        print(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    print(target_pos, query_pos, target_nuc, query_nuc, match, genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                            sequences[gene] = [query_nuc]
                        elif query_pos > results[gene][1]:
                            results[gene][1] = query_pos
                            sequences[gene].append(query_nuc)
                    query_pos += 1
                    if op != 'I':
                        target_pos += 1
            elif op == 'D':
                target_pos += size
                for i in range(size):
                    sequences[gene].append('-')
                print('='*50)
                continue
            elif op == 'I':
                query_pos += size
                print('='*50)
                continue
            print('='*50)
        print('new alignment row'.center(50, '~'))
    return results, sequences

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