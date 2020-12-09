import os
import re
import yaml
import shutil
import subprocess as sp
import pandas as pd
from pathlib import Path
from csv import DictWriter, DictReader
import string
import secrets

from gene_splicer.logger import logger


class Random:

    alphabet = string.ascii_letters

    def __init__(self, size: int = 10) -> None:
        self.size = size
        return

    @staticmethod
    def upper(size=10):
        return ''.join(secrets.choice(self.alphabet)
                       for i in range(size)).upper()


def load_yaml(afile):
    with open(afile) as f:
        return yaml.load(f, Loader=yaml.SafeLoader)


def dump_yaml(data, afile):
    with open(afile, 'w') as o:
        yaml.dump(data, o)


def reverse_and_complement(seq):
    return ''.join(complement_dict[nuc] for nuc in reversed(seq))


def write_annot(adict, filepath):
    """
        adict: {gene (str): [start (int), stop (int)], ...}
    """
    with open(filepath, 'w', newline='') as o:
        fieldnames = ['gene', 'start', 'stop']
        writer = DictWriter(o, fieldnames)
        writer.writeheader()
        for gene, (start, stop) in adict.items():
            writer.writerow({
                'gene': gene.upper(),
                'start': start,
                'stop': stop
            })


def write_fasta(adict, filepath_or_fileobject):
    """Writes a fasta file from a dictionary

    Args:
        adict (dict): A dictionary where the keys are header/sequence names and the values are the sequences
        filepath_or_fileobject (Path/File): A path or file-like object

    Returns:
        Path/File: Path or file-like object
    """
    try:
        with open(filepath_or_fileobject, 'w') as o:
            for key, value in adict.items():
                o.write(f'>{key}\n{value}\n')
    except TypeError:
        for key, value in adict.items():
            filepath_or_fileobject.write(f'>{key}\n{value}\n')
    return filepath_or_fileobject


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
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd'
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
    pattern = re.compile(r'(\d+)([A-Z])')
    cigar = re.findall(pattern, row[5])
    return cigar


class MyContextManager(object):
    def __init__(self, file_name, method):
        self.file_obj = open(file_name, method)

    def __enter__(self):
        return self.file_obj

    def __exit__(self, type, value, traceback):
        self.file_obj.close()


def handle_file_or_fileobject(file_or_fileobject, method):
    try:
        file_or_fileobject.read()
        file_or_fileobject.seek(0)
        return file_or_fileobject
    except AttributeError:
        return MyContextManager(file_or_fileobject, method)


def read_fasta(fasta_file):
    with handle_file_or_fileobject(fasta_file, 'r') as fasta:
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
        if ((annot[gene][0] - offset_start) <= pos <=
            (annot[gene][1] - offset_end)):
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
    newseq = newseq[:pos] + 'G' + newseq[pos + 1:]
    assert newseq[pos] == 'G'

    return newseq


def modify_annot(annot):
    newannot = {}
    for gene, (start, stop) in annot.items():
        if gene in genes_of_interest:
            newannot[gene] = [start, stop]
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
            stop -= 1
        # Convert to 0 base
        newannot[gene] = [start - 1, stop - 1]

    return newannot


def splice_genes(query, target, samfile, annotation):
    results = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        for size, op in row['cigar']:
            size = int(size)
            # logger.debug(f'size: {size}, op: {op}')
            # logger.debug(f'target_pos: {target_pos}, query_pos: {query_pos}')
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                # logger.debug('='*50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                # logger.debug('='*50)
                continue
            elif op in ('M', '=', 'X'):
                for j in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        logger.warning(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    # logger.debug(target_pos, query_pos, target_nuc, query_nuc, match, genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                        elif query_pos > results[gene][1]:
                            results[gene][1] = query_pos
                    query_pos += 1
                    target_pos += 1
            elif op == 'D':
                target_pos += size
                # logger.debug('='*50)
                continue
            elif op == 'I':
                query_pos += size
                # logger.debug('='*50)
                continue
            # logger.debug('='*50)
        # logger.debug('new alignment row'.center(50, '~'))
    return results


def coords_to_genes(coords, query):
    genes = {
        gene: query[coords[0]:coords[1] + 1]
        for gene, coords in coords.items()
    }
    return genes


# This function has not been tested yet
def splice_aligned_genes(query, target, samfile, annotation):
    results = {}
    sequences = {}
    for i, row in samfile.iterrows():
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        for size, op in row['cigar']:
            # print(f'size: {size}, op: {op}')
            # print(f'target_pos: {target_pos}, query_pos: {query_pos}')
            size = int(size)
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                # print('=' * 50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                # print('=' * 50)
                continue
            elif op in ('M', '=', 'X'):
                for i in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        logger.warning(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    # print(target_pos, query_pos, target_nuc, query_nuc, match,
                    #   genes)
                    for gene in genes:
                        if gene not in results:
                            results[gene] = [query_pos, query_pos]
                            sequences[gene] = [query_nuc]
                        elif query_pos > results[gene][1]:
                            results[gene][1] = query_pos
                            sequences[gene].append(query_nuc)
                    query_pos += 1
                    target_pos += 1
            elif op == 'D':
                target_pos += size
                for i in range(size):
                    sequences[gene].append('-')
                # print('=' * 50)
                continue
            elif op == 'I':
                query_pos += size
                # print('=' * 50)
                continue
            # print('=' * 50)
        # print('new alignment row'.center(50, '~'))
    return results, sequences


def load_samfile(samfile_path):
    result = pd.read_table(samfile_path, skiprows=2, header=None)
    result['cigar'] = result.apply(split_cigar, axis=1)
    return result


def minimap2_available(aligner_path='minimap2'):
    cmd = [aligner_path]
    process = sp.run(cmd)
    if process.returncode == 0:
        return True
    return False


# Removes all files in a directory
def clean_dir(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def align(target_seq,
          query_seq,
          outdir=Path(os.getcwd()).resolve(),
          aligner_path='minimap2',
          skip=False):
    outdir = outdir / 'minimap2_aln'
    alignment_path = outdir / 'alignment.sam'
    if skip:
        return alignment_path
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    # Write the query fasta
    query_fasta_path = write_fasta({'query': query_seq},
                                   outdir / 'query.fasta')
    # Write the target fasta
    target_fasta_path = write_fasta({'target': target_seq},
                                    outdir / 'target.fasta')
    cmd = [aligner_path, '-a', target_fasta_path, query_fasta_path]
    with open(alignment_path, 'w') as alignment:
        process = sp.run(cmd, stdout=alignment, errors=True)
    if process.returncode != 0:
        logger.error('Alignment failed! Error: %s' % process.stderr)
        return False
    else:
        return alignment_path


def generate_table_precursor_2(hivseqinr_resultsfile, filtered_file,
                               genes_file, table_precursorfile):
    try:
        seqinr = pd.read_csv(hivseqinr_resultsfile)
    except FileNotFoundError:
        logger.error('hivseqinr could not produce results!')
        # Create an empty results file
        table_precursorfile.touch()
        return False
    # Assign new columns based on split
    # Make sure this matches the join in primer_finder run()
    seqinr[['reference', 'seqtype']] = seqinr['SEQID'].str.split('::',
                                                                 expand=True)

    # Load filtered sequences
    filtered = pd.read_csv(filtered_file)

    # Merge
    merged = seqinr.merge(filtered,
                          left_index=True,
                          right_index=True,
                          how='outer')
    for gene in genes_of_interest:
        merged[gene] = None

    genes_fasta = read_fasta(genes_file)
    genes = dict([x for x in genes_fasta])
    for gene in genes_of_interest:
        try:
            seq = genes[f'>{gene}']
        except KeyError:
            seq = None
        merged[gene] = seq

    # Output csv
    merged[['sequence', 'MyVerdict'] + genes_of_interest].to_csv(
        table_precursorfile, index=False)
    return table_precursorfile


def generate_table_precursor(outpath, table_precursor_path, genes_fasta_path):
    # Load hivseqinr data
    seqinr_path = outpath / 'hivseqinr' / 'Results_Final' / 'Output_MyBigSummary_DF_FINAL.csv'
    try:
        seqinr = pd.read_csv(seqinr_path)
    except FileNotFoundError:
        logger.error('hivseqinr could not produce results!')
        table_precursor_path.touch()
        return False
    # Delete the hxb2 row
    seqinr = seqinr[seqinr['SEQID'] != 'hxb2']
    # Assign new columns based on split
    # Make sure this matches the join in primer_finder run()
    seqinr[['reference', 'seqtype']] = seqinr['SEQID'].str.split('::',
                                                                 expand=True)

    # Load filtered sequences
    filtered_path = outpath / 'filtered.csv'
    filtered = pd.read_csv(filtered_path)

    # Merge
    merged = seqinr.merge(filtered,
                          left_index=True,
                          right_index=True,
                          how='outer')
    for gene in genes_of_interest:
        merged[gene] = None

    data = {}
    genes_fasta = read_fasta(genes_fasta_path)
    genes = dict([x for x in genes_fasta])
    for gene in genes_of_interest:
        try:
            seq = genes[f'>{gene}']
        except KeyError:
            seq = None
        # data.setdefault(gene, []).append(seq)
        data[gene] = seq
    for gene, seq in data.items():
        merged.at[0, gene] = seq

    # Output csv
    merged[['sequence', 'MyVerdict'] + genes_of_interest].to_csv(
        table_precursor_path, index=False)
    return table_precursor_path


def get_softclipped_region(query, alignment):
    size, op = alignment.iloc[0]['cigar'][0]
    if op != 'S':
        logger.warning('Alignment does not start with softclip')
        return
    size = int(size)
    return query[:size]


def sequence_to_coords(query, target, alignment_path, annot):
    aln = load_samfile(alignment_path)
    softclip = get_softclipped_region(query, aln)
    if softclip is None:
        return
    import gene_splicer.probe_finder as probe_finder
    finder = probe_finder.ProbeFinder(softclip, target)
    # query_match = target[finder.start:len(finder.contig_match)]
    target_match = finder.contig_match
    matchlen = len(target_match)
    coords = {}
    for pos in range(finder.start, matchlen - 1):
        genes = get_genes(annot, pos)
        for gene in genes:
            coords.setdefault(gene, [pos, pos])[1] += 1
    return coords


def merge_coords(coords1, coords2):
    """Return new coordinates based on merging coords1 and coords2

    Args:
        coords1 (dict): A dictionary of coords with key=gene_name, value=[start, end]
        coords2 (dict): Same as coords1
    """
    new_coords = {}
    for gene in coords2:
        if gene not in coords1:
            coords1[gene] = coords2[gene][:]
        new_coords[gene] = [
            min(coords2[gene][0], coords1[gene][0]),
            max(coords1[gene][1], coords2[gene][1])
        ]
    return new_coords


## Define some variables
cwd = Path(os.path.realpath(__file__)).parent

# Define genes of interest
genes_of_interest = load_yaml(cwd / 'genes_of_interest.yaml')

# Define some global variables
mixture_dict = {
    'W': 'AT',
    'R': 'AG',
    'K': 'GT',
    'Y': 'CT',
    'S': 'CG',
    'M': 'AC',
    'V': 'AGC',
    'H': 'ATC',
    'D': 'ATG',
    'B': 'TGC',
    'N': 'ATGC',
    '-': 'ATGC'
}

complement_dict = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'W': 'S',
    'R': 'Y',
    'K': 'M',
    'Y': 'R',
    'S': 'W',
    'M': 'K',
    'B': 'V',
    'D': 'H',
    'H': 'D',
    'V': 'B',
    '*': '*',
    'N': 'N',
    '-': '-'
}

hxb2 = next(read_fasta(cwd / 'hxb2.fasta'))[1]
mod_hxb2 = modify_reference(hxb2)

annot = {
    x['gene']: [int(x['start']), int(x['stop'])]
    for x in read_csv(cwd / 'annot.csv')
}
mod_annot = modify_annot(annot)