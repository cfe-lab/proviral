import re
import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
import os
import pandas as pd
import subprocess
import requests
import io
import zipfile
import shutil
from pathlib import Path
import sys
import logging

from gene_splicer.logger import logger
from gene_splicer.primer_finder_class import PrimerFinder
from gene_splicer.outcome_summary import OutcomeSummary
from gene_splicer.hivseqinr import Hivseqinr
import gene_splicer.utils as utils

mixture_dict = utils.mixture_dict
reverse_and_complement = utils.reverse_and_complement

logger = logging.getLogger('gene_splicer')


# Note these are 1-based indicies
primers = {
    'fwd': {
        'seq': 'GCGCCCGAACAGGGACYTGAAARCGAAAG',
        'nomix': 'GCGCCCGAACAGGGACCTGAAAGCGAAAG',
        # convert to 0-base index
        'hxb2_start': 638 - 1,
        'hxb2_end': 666
    },
    'rev': {
        'seq': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        'nomix': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
        # convert to 0-base index
        'hxb2_start': 9604 - 1,
        'hxb2_end': 9632
    }
}


def parse_args():
    parser = ArgumentParser(
        description='Search sequences for primers and identify errors',
        formatter_class=ArgumentDefaultsHelpFormatter)

    # Inputs
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences',
                        type=FileType())
    parser.add_argument('conseqs_csv',
                        help='CSV file with conseq sequences',
                        type=FileType())
    parser.add_argument(
        '-s',
        '--sample_size',
        type=int,
        help=
        'Length of sequence (probe) from each end of sample to search for primer',
        default=50)
    parser.add_argument(
        '-e',
        '--extended_size',
        type=int,
        help=
        'Length of sequence to look for the entire primer if only partial is found',
        default=200)

    # Outputs
    parser.add_argument(
        'table_precursor_csv',
        help=
        'CSV file containing gene regions and hivseqinr verdict for sample',
        type=Path)

    parser.add_argument(
        'aligned_table_precursor_csv',
        help=
        'CSV file containing aligned gene regions to hxb2 and hivseqinr verdict for sample',
        type=Path)

    parser.add_argument(
        'outcome_summary_csv',
        help=
        'CSV file containing debugging information for samples',
        type=Path)

    parser.add_argument('-o',
                        '--outpath',
                        help='The path to save the output',
                        default=Path(os.getcwd()),
                        type=Path)
    parser.add_argument('--disable_hivseqinr',
                        action='store_true',
                        help='Disable running hivseqinr')
    parser.add_argument(
        '--nodups',
        action='store_false',
        help=
        'Set this flag to disable the removal of duplicate samples that pass QC'
    )
    parser.add_argument(
        '--skip_align',
        action='store_true',
        help=
        'For testing only. Enable to disable alignment calls (minimap does not run on windows)'
    )
    return parser.parse_args()


def make_path(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def find_primers(csv_filepath,
                 outpath,
                 seqtype,
                 sample_size=50,
                 extended_size=200):
    make_path(outpath)
    v3_reference = 'HIV1-CON-XX-Consensus-seed'
    columns = [
        'reference', 'error', 'sequence', 'seqlen', 'nmixtures', 'seqtype'
    ]
    for direction in primers:
        for column_type in [
                'error', 'canonical_primer_seq', 'sample_primer_seq',
                'sample_primer_start', 'sample_primer_end',
                'sample_primer_size', 'hxb2_sample_primer_start',
                'hxb2_sample_primer_end'
        ]:
            columns.append(direction + '_' + column_type)
    non_tcga = re.compile(r'[^TCGA-]+')
    outfilepath = outpath / f'{seqtype}_primer_analysis.csv'
    outfile = open(outfilepath, 'w')
    writer = DictWriter(outfile, columns, lineterminator='\n')
    writer.writeheader()
    reader = DictReader(csv_filepath)
    # projects = ProjectConfig.loadDefault()
    # hxb2 = projects.getReference('HIV1-B-FR-K03455-seed')
    hxb2 = utils.hxb2
    skipped = {}
    total = 0
    viable = 0
    contig_num = 0

    for row in reader:
        total += 1
        seed_name = row.get('genotype') or row.get('ref') or row['region']
        conseq_cutoff = row.get('consensus-percent-cutoff')
        contig_num += 1
        contig_name = f'{contig_num}_{seed_name}'
        uname = f'{contig_name}_{contig_num}'

        new_row = dict(reference=contig_name, seqtype=seqtype)
        contig_seq: str = row.get('contig') or row['sequence']
        contig_seq = contig_seq.upper()
        new_row['seqlen'] = len(contig_seq)
        new_row['sequence'] = contig_seq
        # print(uname)
        if conseq_cutoff and conseq_cutoff != 'MAX':
            skipped[uname] = 'contig not MAX'
            new_row['error'] = skipped[uname]
            writer.writerow(new_row)
            continue

        # I don't think we need this try block if I am already removing non-proviral stuff up on line 218
        try:
            # If "region" is a column of row, then we are looking at a conseq and not a contig. Only conseqs can have V3 sequences so if we can't access this key we do nothing
            if row['region'] == v3_reference:
                skipped[uname] = 'is V3 sequence'
                new_row['error'] = skipped[uname]
                writer.writerow(new_row)
                continue
        except KeyError:
            pass

        # Determine if sequence has internal Xs
        x_locations = [i for i, j in enumerate(contig_seq) if j == 'X']
        if any([(sample_size < i < len(contig_seq) - (sample_size))
                for i in x_locations]):
            skipped[uname] = 'contig sequence contained internal X'
            new_row['error'] = skipped[uname]
            writer.writerow(new_row)
            continue
        found_non_tcga = re.findall(non_tcga, contig_seq)
        mixtures = len([x for x in found_non_tcga if x[0].upper() != 'X'])
        if (mixtures > 1):
            skipped[uname] = 'contig sequence contained non-TCGA/gap'
            new_row['error'] = skipped[uname]
            new_row['nmixtures'] = mixtures
            writer.writerow(new_row)
            continue
        prime5_seq = contig_seq[:sample_size]
        prime3_seq = contig_seq[-sample_size:]
        for key in columns:
            if key not in [
                    'contig', 'seqlen', 'error', 'sequence', 'reference',
                    'seqtype'
            ]:
                new_row[key] = None
        for end, seq in ((5, prime5_seq), (3, prime3_seq)):
            if end == 5:
                direction = 'fwd'
                hxb2_target_start = primers[direction]['hxb2_start']
                hxb2_target_end = primers[direction]['hxb2_end'] + 100
            else:
                direction = 'rev'
                hxb2_target_start = primers[direction]['hxb2_start'] - 100
                hxb2_target_end = primers[direction]['hxb2_end']

            prefix = direction + '_'

            if 'X' in seq:
                seqlen = len(seq)
                oldseq = seq
                seq = handle_x(seq)
                if not seq or len(seq) < seqlen / 6:
                    new_row[prefix + 'error'] = 'too many X in sequence'
                    continue

            # if 'X' in seq:
            #     import pdb; pdb.set_trace()
            finder = PrimerFinder(contig_seq,
                                  primers[direction]['seq'],
                                  direction,
                                  hxb2_target_start,
                                  hxb2_target_end,
                                  sample_size=sample_size)

            if not finder.is_full_length:
                finder2 = PrimerFinder(contig_seq,
                                       primers[direction]['seq'],
                                       direction,
                                       hxb2_target_start,
                                       hxb2_target_end,
                                       sample_size=extended_size)
                if finder2.is_full_length:
                    finder = finder2

            # Natalie's request
            # If a primer is not found at all, have a custom error for it, if there is something found but it did not pass secondary validation then make a different error for that
            if not finder.sample_primer:
                new_row[prefix + 'error'] = 'primer was not found'
                continue
            # This can happen if there are too many combinations of the sequence to unpack (too many mixtures or dashes etc)
            if not finder.is_valid:
                new_row[prefix + 'error'] = 'primer was not found'
                continue
            new_row[prefix +
                    'canonical_primer_seq'] = primers[direction]['seq']
            new_row[prefix + 'sample_primer_seq'] = finder.sample_primer
            new_row[prefix + 'sample_primer_start'] = finder.start
            new_row[prefix + 'sample_primer_end'] = finder.end
            new_row[prefix + 'sample_primer_size'] = len(finder.sample_primer)
            new_row[prefix + 'hxb2_sample_primer_start'] = finder.hxb2_start
            new_row[prefix + 'hxb2_sample_primer_end'] = finder.hxb2_end
        writer.writerow(new_row)
        viable += 1

    # If no reads remapped, contig/conseq does not exist, write it as an error
    if total == 0:
        new_row = dict(reference=None,
                       error='No contig/conseq constructed',
                       sequence=None,
                       seqlen=None,
                       nmixtures=None)
        for prefix in ('fwd_', 'rev_'):
            new_row[prefix + 'error'] = None
            new_row[prefix + 'canonical_primer_seq'] = None
            new_row[prefix + 'sample_primer_seq'] = None
            new_row[prefix + 'sample_primer_start'] = None
            new_row[prefix + 'sample_primer_end'] = None
            new_row[prefix + 'sample_primer_size'] = None
            new_row[prefix + 'hxb2_sample_primer_start'] = None
            new_row[prefix + 'hxb2_sample_primer_end'] = None
        writer.writerow(new_row)

    outfile.close()
    return outfilepath


def handle_x(sequence):
    length = len(sequence)
    midpoint = length / 2
    x_positions = [i for i, j in enumerate(sequence) if j == 'X']
    try:
        rightmost_x = max([x for x in x_positions if x <= midpoint])
    except ValueError:
        rightmost_x = None
    try:
        leftmost_x = min([x for x in x_positions if x > midpoint])
    except ValueError:
        leftmost_x = None
    if rightmost_x is not None and leftmost_x is not None:
        sequence = sequence[rightmost_x + 1:leftmost_x]
    elif rightmost_x is not None:
        sequence = sequence[rightmost_x + 1:]
    elif leftmost_x is not None:
        sequence = sequence[:leftmost_x]
    return sequence


def validate_primer(finder, finder_seq, target, hxb2_target, tolerance=1):
    if finder.is_reversed:
        finder_seq = reverse_and_complement(finder_seq)
    error = None
    primer_in_finder = None
    finder_seqsize = len(finder_seq)
    matched_finder_size = len(finder.contig_match)
    overhang = matched_finder_size - finder_seqsize
    primer_in_finder_hxb2_start_coord = max(finder.start, target['hxb2_start'])
    primer_in_finder_hxb2_end_coord = min(finder.start + finder_seqsize,
                                          target['hxb2_end'])
    # Get the coordinates of the primer relative to the finder sequence
    primer_in_finder_start_coord = primer_in_finder_hxb2_start_coord - finder.start
    primer_in_finder_end_coord = primer_in_finder_hxb2_end_coord - finder.start
    if ((target['hxb2_start'] == 9603) and
        (primer_in_finder_end_coord >= matched_finder_size - overhang)):
        primer_in_finder_start_coord -= overhang
        primer_in_finder_end_coord -= overhang

    # Get the primer sequence of the finder sequence
    primer_in_finder = finder_seq[
        primer_in_finder_start_coord:primer_in_finder_end_coord]

    # Get the sequence of the true primer that overlaps the finder sequence
    primer_start_coord = max(
        0, primer_in_finder_hxb2_start_coord - target['hxb2_start'])
    primer_end_coord = primer_start_coord + len(primer_in_finder)
    real_primer = target['seq'][primer_start_coord:primer_end_coord]
    result = {
        'finder_seq': finder_seq,
        'target_seq': primer_in_finder,
        'contig_hxb2_start': finder.start,
        'contig_hxb2_end': matched_finder_size,
        'seq_start': primer_in_finder_start_coord,
        'seq_end': primer_in_finder_end_coord,
        'hxb2_start': primer_in_finder_hxb2_start_coord,
        'hxb2_end': primer_in_finder_hxb2_end_coord,
        'overhang': overhang,
        'contig_match': finder.contig_match,
        'dist': 0,
        'finder_dist': finder.dist,
        'real_primer': real_primer,
        'error': error
    }

    if len(real_primer) != len(primer_in_finder):
        result[
            'error'] = 'Real primer did not match length of primer in finder'
        return result
    if not real_primer:
        result['error'] = 'real primer not found at expected coordinates'
        return result
    elif not primer_in_finder:
        result[
            'error'] = 'primer in contig sequence not found at expected coordinates'
        return result
    for i in range(len(real_primer)):
        try:
            our_nuc = primer_in_finder[i]
        except IndexError:
            pass
        real_nuc = real_primer[i]
        if our_nuc != real_nuc:
            if real_nuc in mixture_dict and our_nuc in mixture_dict[real_nuc]:
                pass
            else:
                result['dist'] += 1
        if result['dist'] > tolerance:
            result['error'] = 'mismatches in primer > tolerance'
            return result
    # if finder.dist > (finder_seqsize / 10):
    #     result['error'] = f'Poor alignment (finder_dist > {finder_seqsize/10})'
    #     return result
    return result


def load_csv(csv_filepath, seqtype, results=None):
    if results is None:
        results = {}
    df = pd.read_csv(csv_filepath)
    df['seqtype'] = seqtype
    results[seqtype] = df
    return results


def add_primers(row):
    # Strip the primers out
    newseq = row.sequence[int(row.fwd_in_probe_size
                              ):-int(row.rev_in_probe_size)]
    # Add the primers in
    #     newseq = primers['fwd']['seq'] + newseq + primers['rev']['seq']
    row.sequence = newseq
    return row


def remove_primers(row):
    # Strip the primers out
    newseq = row.sequence[int(row.fwd_sample_primer_size
                              ):-int(row.rev_sample_primer_size)]
    row.sequence = newseq
    return row


def filter_df(df, nodups=True):
    filtered = df[(df['error'].isna()
                   & df['fwd_error'].isna()
                   & df['rev_error'].isna())]
    filtered = filtered.apply(remove_primers, axis=1)
    if nodups:
        if len(filtered) > 1:
            # Return empty dataframe
            return pd.DataFrame(columns=df.columns, index=df.index)
    # Remove any rows with references containing "reverse" or "unknown"
    filtered = filtered[(~filtered['reference'].str.contains('reverse'))
                        & (~filtered['reference'].str.contains('unknown'))]
    filtered = filtered[['reference', 'sequence', 'seqtype']]
    return filtered


def run(contigs_csv,
        conseqs_csv,
        outpath,
        outcome_summary_path,
        disable_hivseqinr=False,
        nodups=True,
        sample_size=50,
        extended_size=200):
    contigs_out = find_primers(contigs_csv,
                               outpath,
                               'contigs',
                               sample_size=sample_size,
                               extended_size=extended_size)
    conseqs_out = find_primers(conseqs_csv,
                               outpath,
                               'conseqs',
                               sample_size=sample_size,
                               extended_size=extended_size)
    dfs = load_csv(contigs_out, 'contigs')
    dfs = load_csv(conseqs_out, 'conseqs', dfs)
    files = []
    contigs_df = dfs['contigs'].fillna('')
    conseqs_df = dfs['conseqs'].fillna('')
    # Generate outcome summary
    OutcomeSummary(contigs_df, conseqs_df, outcome_summary_path)
    filtered_contigs = filter_df(contigs_df, nodups)
    filtered_conseqs = filter_df(conseqs_df, nodups)
    joined = filtered_contigs.merge(filtered_conseqs,
                                    left_index=True,
                                    right_index=True,
                                    suffixes=('_contig', '_conseq'),
                                    how='outer')
    joined['sequence'] = joined['sequence_conseq'].fillna(
        joined['sequence_contig'])
    joined['seqtype'] = joined['seqtype_conseq'].fillna(
        joined['seqtype_contig'])
    joined['seqlen'] = joined['sequence'].fillna('').str.len()
    joined['reference'] = joined['reference_conseq'].fillna(
        joined['reference_contig'])
    joined = joined[['reference', 'seqtype', 'sequence', 'seqlen']]
    if joined['sequence'].isnull().all():
        logger.warning('No contigs and no conseqs passed QC, exiting')
        return files
    joined.to_csv(outpath / 'filtered.csv', index=False)
    synthetic_primers_added_fasta_path = outpath / 'synthetic_primers.fasta'
    no_primers_fasta_path = outpath / 'no_primers.fasta'
    if joined['sequence'].isnull().all():
        if not disable_hivseqinr:
            hivseqinr = Hivseqinr(outpath / 'hivseqinr',
                                  synthetic_primers_added_fasta_path)
        files.append(no_primers_fasta_path)
        return files
    with open(synthetic_primers_added_fasta_path,
              'w') as o, open(no_primers_fasta_path, 'w') as o2:
        # Write hxb2 as positive control
        o.write(f'>hxb2\n{utils.hxb2}\n')
        o2.write(f'>hxb2\n{utils.hxb2}\n')
        for row in joined.itertuples():
            # I don't remember why it was necessary to replace dashes with underscores but I think it was because HIVSEQINR doesn't like dashes in names
            # I've commented it out for now
            # header = f'>{row.name}_{row.sample}_{row.reference}_{row.seqtype}'.replace('-', '_')
            # The header delimiter, this must match the split in gene_splicer
            header = '>' + '::'.join((row.reference, row.seqtype))
            o.write(
                f'{header}\n{primers["fwd"]["nomix"] + row.sequence.replace("-", "") + primers["rev"]["nomix"]}\n'
            )
            o2.write(f'{header}\n{row.sequence.replace("-", "")}\n')
    if not disable_hivseqinr:
        hivseqinr = Hivseqinr(outpath / 'hivseqinr',
                              synthetic_primers_added_fasta_path)
    files.append(no_primers_fasta_path)
    return files


def main():
    args = parse_args()
    fasta_files = run(contigs_csv=args.contigs_csv,
                      conseqs_csv=args.conseqs_csv,
                      outpath=args.outpath.resolve(),
                      outcome_summary_path=args.outcome_summary_csv,
                      disable_hivseqinr=args.disable_hivseqinr,
                      nodups=args.nodups,
                      sample_size=args.sample_size,
                      extended_size=args.extended_size)
    return {'fasta_files': fasta_files, 'args': args}


if __name__ in ('__main__', '__live_coding__'):
    main()
