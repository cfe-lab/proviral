import csv
import logging
import os
import re
import typing
from typing import TextIO, Mapping, Dict, Set, List, Iterable, Tuple, Optional, Literal

import yaml
import json
import shutil
from cfeproviral.version import get_version
import subprocess as sp
import pandas as pd
import glob
from pathlib import Path
from csv import DictWriter, DictReader
from itertools import groupby
from operator import itemgetter

logger = logging.getLogger(__name__)

Backend = Optional[Literal["CFEIntact", "HIVSeqinR"]]

LEFT_PRIMER_END = 666
RIGHT_PRIMER_START = 9604


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


def split_cigar(string):
    pattern = re.compile(r'(\d+)([A-Z])')
    cigar = re.findall(pattern, string)
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
    # newannot = dict(annot)
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
    for i, row in enumerate(samfile):
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        cigar = row[5]
        for size, op in split_cigar(cigar):
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


def splice_aligned_genes(query, target, samfile, annotation):
    results = {}
    sequences = {}
    for i, row in enumerate(samfile):
        # Subtract 1 to convert target position to zero-base
        target_pos = int(row[3]) - 1
        query_pos = None
        cigar = row[5]
        for size, op in split_cigar(cigar):
            # print(f'size: {size}, op: {op}')
            # print(f'target_pos: {target_pos}, query_pos: {query_pos}')
            size = int(size)
            # If the first section is hard-clipped the query should start at
            # the first non-hard-clipped-based. The target should also be offset
            # by the size of the hard-clip
            if op == 'H' and query_pos is None:
                query_pos = size
                print('=' * 50)
                continue
            elif query_pos is None:
                query_pos = 0
            if op == 'S':
                query_pos += size
                print('=' * 50)
                continue
            elif op in ('M', '=', 'X'):
                for _ in range(size):
                    try:
                        target_nuc = target[target_pos]
                    except IndexError:
                        logger.warning(f'{target_pos} not in range of target')
                        break
                    query_nuc = query[query_pos]
                    match = (target_nuc == query_nuc)
                    genes = get_genes(annotation, target_pos)
                    print(target_pos, query_pos, target_nuc, query_nuc, match,
                          genes)
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
                for _ in range(size):
                    sequences[gene].append('-')
                print('=' * 50)
                continue
            elif op == 'I':
                query_pos += size
                print('=' * 50)
                continue
            print('=' * 50)
        print('new alignment row'.center(50, '~'))
    return results, sequences


def load_samfile(samfile_path: Path) -> List[List[str]]:
    with open(samfile_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')

        result = []
        for row in reader:
            # Skip header lines (lines starting with '@')
            if row[0].startswith('@'):
                continue
            result.append(row)

    return result


def aligner_available(aligner_path='minimap2'):
    cmd = [aligner_path, '-h']
    try:
        process = sp.run(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
        if process.returncode == 0:
            return True
    except Exception:
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
          query_name,
          outdir=Path(os.getcwd()).resolve(),
          aligner_path='minimap2'):
    if not aligner_available(aligner_path):
        raise FileNotFoundError(f'Aligner {aligner_path} not found.')
    outdir = outdir / query_name
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    # Write the query fasta
    query_fasta_path = write_fasta({query_name: query_seq},
                                   outdir / 'query.fasta')
    # Write the target fasta
    target_fasta_path = write_fasta({'MOD_HXB2': target_seq},
                                    outdir / 'target.fasta')
    cmd = [aligner_path, '-a', target_fasta_path, query_fasta_path]
    alignment_path = outdir / 'alignment.sam'
    log_path = outdir / 'minimap2.log'
    with alignment_path.open('w') as alignment, log_path.open('w') as log_file:
        process = sp.run(cmd,
                         stdout=alignment,
                         stderr=log_file,
                         check=True)
    if process.returncode != 0:
        logger.error('Alignment failed! Details in %s.', log_path)
        return False
    else:
        return alignment_path

CFEINTACT_ERRORS_TABLE = [
    'UnknownNucleotide',
    'NonHIV',
    'LongDeletion',
    'InternalInversion',
    'Scramble',
    'APOBECHypermutation',
    'MajorSpliceDonorSiteMutated',
    'PackagingSignalDeletion',
    'PackagingSignalNotComplete',
    'RevResponseElementDeletion',
    'MisplacedORF',
    'WrongORFNumber',
    'Deletion',
    'Insertion',
    'InternalStop',
    'Frameshift',
    'MutatedStopCodon',
    'MutatedStartCodon',
    'SequenceDivergence',
    ]

def iterate_cfeintact_verdicts_1(directory: Path, intact: Set[str] = set()) -> Iterable[Tuple[str, str]]:
    intact = set()

    def get_verdict(SEQID: str, all_defects) -> Tuple[str, str]:
        if all_defects:
            ordered = sorted(all_defects, key=CFEINTACT_ERRORS_TABLE.index)
            verdict = ordered[0]
        else:
            verdict = "Intact"

        return (SEQID, verdict)

    with open(os.path.join(directory, 'holistic.csv'), 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["intact"] == "True":
                intact.add(row["qseqid"])
                SEQID = row["qseqid"]
                yield get_verdict(SEQID, all_defects=[])

    with open(os.path.join(directory, 'defects.csv'), 'r') as f:
        reader = csv.DictReader(f)
        grouped = groupby(reader, key=itemgetter('qseqid'))
        for sequence_name, defects in grouped:
            if sequence_name not in intact:
                all_defects = [defect['code'] for defect in defects]
                yield get_verdict(sequence_name, all_defects)


def iterate_cfeintact_verdicts(outpath: Path) -> Iterable[Tuple[str, str]]:
    intact: Set[str] = set()

    for directory in outpath.glob('cfeintact*'):
        yield from iterate_cfeintact_verdicts_1(directory, intact)


def get_cfeintact_verdicts(name, outpath):
    column_names = ['SEQID', 'MyVerdict']
    data = iterate_cfeintact_verdicts(outpath)
    return pd.DataFrame(data, columns=column_names)


def iterate_hivseqinr_verdicts_1(directory: Path) -> Iterable[Tuple[str, str]]:
    path = directory / "Results_Final" / "Output_MyBigSummary_DF_FINAL.csv"
    if not path.is_file():
        logger.error("Missing HIVSeqinR result in %r.", str(path))
        return

    with path.open() as fd:
        reader = csv.DictReader(fd)
        for row in reader:
            yield (row["SEQID"], row["MyVerdict"])


def iterate_hivseqinr_verdicts(outpath: Path) -> Iterable[Tuple[str, str]]:
    for directory in outpath.glob('hivseqinr*'):
        yield from iterate_hivseqinr_verdicts_1(directory)


def get_hivseqinr_verdicts(name, outpath):
    column_names = ['SEQID', 'MyVerdict']
    data = iterate_hivseqinr_verdicts(outpath)
    return pd.DataFrame(data, columns=column_names)


def generate_table_precursor(name, outpath, add_columns=None):
    # Output csv
    precursor_path: Path = outpath / 'table_precursor.csv'

    # Load filtered sequences
    filtered_path = outpath / (name + '_filtered.csv')
    filtered = pd.read_csv(filtered_path)
    # Load hivseqinr data or Cfeintact results

    if any(outpath.glob('cfeintact*')):
        results = get_cfeintact_verdicts(name, outpath)
    elif any(outpath.glob('hivseqinr*')):
        results = get_hivseqinr_verdicts(name, outpath)
    else:
        raise RuntimeError("Neither CFEIntact nor HIVSeqinR directory exists.")

    try:
        # Assign new columns based on split
        results[['name', 'sample', 'reference',
                'seqtype']] = results['SEQID'].str.split('::', expand=True)
        # Merge
        merged = results.merge(filtered, on='sample')
    except ValueError:
        with precursor_path.open('w') as output_file:
            writer = DictWriter(output_file,
                                ['sample', 'sequence', 'MyVerdict'] +
                                genes_of_interest)
            writer.writeheader()
        return precursor_path

    for gene in genes_of_interest:
        merged[gene] = None

    data = {}
    for index, row in merged.iterrows():
        folder = outpath / row['sample']
        genes_fasta_path = folder / 'genes.fasta'
        if not os.path.isfile(genes_fasta_path):
            print(f'No genes for {row["sample"]}')
            continue
        genes_fasta = read_fasta(genes_fasta_path)
        genes = dict([x for x in genes_fasta])
        for gene in genes_of_interest:
            try:
                seq = genes[f'>{gene}']
            except KeyError:
                seq = None
            data.setdefault(gene, []).append(seq)
    for gene, seqs in data.items():
        merged[gene] = seqs

    if add_columns:
        for key, val in add_columns.items():
            merged[key] = val
    merged['cfeproviral_version'] = get_version()
    if not results.empty:
        merged[['sample', 'sequence', 'MyVerdict'] + genes_of_interest + ['cfeproviral_version']].to_csv(
            precursor_path, index=False)
    else:
        merged[['sample', 'sequence'] + genes_of_interest + ['cfeproviral_version']].to_csv(precursor_path,
                                                                  index=False)
    return precursor_path


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
    merged['cfeproviral_version'] = get_version()
    merged[['sequence', 'MyVerdict'] + genes_of_interest + ['cfeproviral_version']].to_csv(
        table_precursorfile, index=False)
    return table_precursorfile


def get_softclipped_region(query, alignment, alignment_path):
    try:
        first_match = alignment[0]
    except IndexError:
        logger.warning('No alignment in %s!', alignment_path)
        return

    cigar = first_match[5]
    if cigar == '*':
        logger.warning('No alignment in %s!', alignment_path)
        return

    size, op = split_cigar(cigar)[0]
    if op != 'S':
        return

    size = int(size)
    return query[:size]


def sequence_to_coords(query, target, alignment_path, annotations):
    aln = load_samfile(alignment_path)
    softclip = get_softclipped_region(query, aln, alignment_path)
    if softclip is None:
        return
    import cfeproviral.probe_finder
    finder = cfeproviral.probe_finder.ProbeFinder(softclip, target)
    if not finder.valid:
        return None
    # query_match = target[finder.start:len(finder.contig_match)]
    target_match = finder.contig_match
    matchlen = len(target_match)
    coords = {}
    for pos in range(finder.start, matchlen - 1):
        genes = get_genes(annotations, pos)
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
            new_coords[gene] = coords2[gene][:]
            continue
        new_coords[gene] = [
            min(coords2[gene][0], coords1[gene][0]),
            max(coords1[gene][1], coords2[gene][1])
        ]
    for gene in coords1:
        if gene not in coords2:
            new_coords[gene] = coords1[gene][:]
            continue
        new_coords[gene] = [
            min(coords1[gene][0], coords2[gene][0]),
            max(coords2[gene][1], coords1[gene][1])
        ]
    return new_coords


def filter_valid(df):
    # Remove any row that has no errors
    filtered = df[(~df['error'].isna())
                  | (~df['fwd_error'].isna())
                  | (~df['rev_error'].isna())].copy()
    # Remove contig not max errors and V3 errors
    filtered = filtered[(filtered['error'] != 'contig not MAX')
                        & (filtered['error'] != 'is V3 sequence')]
    # Set error field for duplicates
    #filtered.loc[filtered.duplicated(subset='sample', keep=False),
    #             'error'] = 'duplicate'
    filtered = filtered[(~filtered['reference'].str.contains('reverse'))
                        & (~filtered['reference'].str.contains('unknown'))]
    return filtered


def genFailureSummary(contigs_df, conseqs_df, outpath):
    filtered_contigs = filter_valid(contigs_df)
    filtered_conseqs = filter_valid(conseqs_df)
    contigs_simple = filtered_contigs[[
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error'
    ]]
    conseqs_simple = filtered_conseqs[[
        'sample', 'run_name', 'reference', 'error', 'fwd_error', 'rev_error'
    ]]
    concat = pd.concat([contigs_simple, conseqs_simple])
    outfile = outpath / 'failure_summary.csv'
    concat.to_csv(outfile, index=False)
    return outfile


# Checks if a value is nan, handles strings
def isNan(num):
    return num != num


def get_samples_from_cascade(cascade_csv: typing.IO,
                             default_sample_name: Optional[str] = None):
    all_samples = {}
    reader = DictReader(cascade_csv)
    if reader.fieldnames and 'sample' in reader.fieldnames:
        for row in reader:
            all_samples[row['sample']] = {
                'remap': int(row['remap']),
                'micall_version': row.get('micall_version')
            }
        return all_samples
    rows = list(reader)
    assert len(rows) == 1, len(rows)
    remap_count = int(rows[0]['remap'])
    micall_version = rows[0].get('micall_version')
    return {default_sample_name: {'remap': remap_count, 'micall_version': micall_version}}


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

# Note these are 1-based indicies
primers = {
    'fwd': {
        'seq': 'GCGCCCGAACAGGGACYTGAAARCGAAAG',
        'nomix': 'GCGCCCGAACAGGGACCTGAAAGCGAAAG',
        # convert to 0-base index
        'hxb2_start': 638 - 1,
        'hxb2_end': 666,
        'direction': 'fwd'
    },
    'rev': {
        'seq': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        'nomix': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        # convert to 0-base index
        'hxb2_start': 9604 - 1,
        'hxb2_end': 9632,
        'direction': 'rev'
    }
}
