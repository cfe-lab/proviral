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
import math

import Levenshtein
from gotoh import align_it
from gene_splicer.logger import logger
from gene_splicer.primer_finder_errors import PrimerFinderErrors
from gene_splicer.primer_finder_class import PrimerFinder
from gene_splicer.outcome_summary import OutcomeSummary
from gene_splicer.hivseqinr import Hivseqinr
from gene_splicer.helpers.proviral_helper import ProviralHelper
from gene_splicer.probe_finder import ProbeFinder

import gene_splicer.utils as utils
mixture_dict = utils.mixture_dict
reverse_and_complement = utils.reverse_and_complement

# from micall.core.project_config import ProjectConfig
# from micall.utils.translation import mixture_dict, reverse_and_complement

logger = logging.getLogger('gene_splicer')

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


def parse_args():
    parser = ArgumentParser(
        description='Search sequences for primers and identify errors',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences',
                        type=FileType())
    parser.add_argument('conseqs_csv',
                        help='CSV file with conseq sequences',
                        type=FileType())
    parser.add_argument('cascade_csv',
                        help='CSV file from MiCall output',
                        type=FileType())
    parser.add_argument('name', help='A name for the analysis')
    parser.add_argument(
        '-p',
        '--sample_size',
        type=int,
        help=
        'Length of sequence (probe) from each end of sample to search for primer',
        default=50)
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
        '--split',
        type=int,
        default=1,
        help=
        'To avoid memory issues in hivseqinr, split the resulting qc-passed sequences into this number of fastas, each will be processed sequentially and then all will be merged into the final result'
    )
    return parser.parse_args()


def make_path(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def find_primers(
        csv_filepath,
        outpath,
        run_name,
        all_samples,
        # The sequence type, can be either 'contigs' or 'conseqs'
        seqtype,
        sample_size=50,
        extended_length=200,
        test=False):
    reversed_conseqs = {}
    proviral_helper = ProviralHelper(force_all_proviral=test)
    errors = PrimerFinderErrors()
    v3_reference = 'HIV1-CON-XX-Consensus-seed'
    print(run_name)
    make_path(outpath)
    columns = [
        'run_name', 'sample', 'reference', 'error', 'sequence', 'seqlen',
        'nmixtures'
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
    outfilepath = os.path.join(outpath, f'{run_name}_primer_analysis.csv')
    outfile = open(outfilepath, 'w')
    writer = DictWriter(outfile, columns, lineterminator='\n')
    writer.writeheader()
    reader = DictReader(csv_filepath)

    # all_samples is a dict whose keys are sample names and values are the remap column (constructed from cascade.csv file)
    for sample in all_samples:
        # If no reads remapped, contig/conseq does not exist, write it as an error
        if all_samples[sample] == 0:
            new_row = dict(run_name=run_name,
                           sample=sample,
                           reference=None,
                           error=errors.no_sequence,
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

    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):

        # Do not analyze non-proviral samples
        if not proviral_helper.is_proviral(sample_name):
            logger.debug('Skipping sample "%s" because it is non-proviral' %
                         sample_name)
            continue

        contig_num = 0

        for row in sample_rows:
            contig_num += 1
            seed_name = row.get('ref') or row['region'] or row.get('genotype')

            conseq_cutoff = row.get('consensus-percent-cutoff')
            contig_name = f'{contig_num}-{seed_name}'

            new_row = dict(run_name=run_name,
                           sample=sample_name,
                           reference=contig_name)

            contig_seq: str = row.get('contig') or row['sequence']
            contig_seq = contig_seq.upper()

            # Determine if sequence has internal Xs
            x_locations = [i for i, j in enumerate(contig_seq) if j == 'X']
            if any([(sample_size < i < len(contig_seq) - (sample_size))
                    for i in x_locations]):
                new_row['error'] = errors.low_internal_cov
                writer.writerow(new_row)
                continue

            # Skip v3 sequences
            try:
                # If "region" is a column of row, then we are looking at a conseq and not a contig. Only conseqs can have V3 sequences so if we can't access this key we do nothing
                if row['region'] == v3_reference:
                    new_row['error'] = errors.non_proviral
                    writer.writerow(new_row)
                    continue
            except KeyError:
                pass

            # Keep track of reversed conseq seeds
            if seqtype == 'conseqs':
                conseq_num, seed = seed_name.split('-', 1)
                conseq_num = int(conseq_num)
                if 'reversed' in seed_name:
                    reversed_conseqs[f'{conseq_num}-{sample_name}'] = seed

            # If corresponding conseq is reversed, reverse contig
            # TODO Reverse conseq as well if reversed seed
            if f'{contig_num}-{sample_name}' in reversed_conseqs:
                contig_seq = utils.reverse_and_complement(contig_seq)

            new_row['seqlen'] = len(contig_seq)
            new_row['sequence'] = contig_seq

            # If the sample is non-proviral, skip
            if not proviral_helper.is_proviral(row.get('sample')):
                new_row['error'] = errors.non_proviral
                writer.writerow(new_row)
                continue

            # If percent consensus cutoff is not max, skip
            if conseq_cutoff and conseq_cutoff != 'MAX':
                new_row['error'] = errors.not_max
                writer.writerow(new_row)
                continue

            # Determine if sequence has non-tcga characters
            found_non_tcga = re.findall(non_tcga, contig_seq)
            mixtures = len([x for x in found_non_tcga if x[0].upper() != 'X'])
            if (mixtures > 1):
                new_row['error'] = errors.non_tcga
                new_row['nmixtures'] = mixtures
                writer.writerow(new_row)
                continue

            prime5_seq = contig_seq[:sample_size]
            extended_prime5_seq = contig_seq[:extended_length]
            prime3_seq = contig_seq[-sample_size:]
            extended_prime3_seq = contig_seq[-extended_length:]
            for key in columns:
                if key not in [
                        'sample', 'contig', 'seqlen', 'error', 'sequence',
                        'run_name', 'reference'
                ]:
                    new_row[key] = None

            for end, seq, extended_seq in [
                (5, prime5_seq, extended_prime5_seq),
                (3, prime3_seq, extended_prime3_seq)
            ]:
                if end == 5:
                    direction = 'fwd'
                    hxb2_target_start = primers[direction]['hxb2_start']
                    hxb2_target_end = primers[direction]['hxb2_end']
                else:
                    direction = 'rev'
                    hxb2_target_start = primers[direction]['hxb2_start']
                    hxb2_target_end = primers[direction]['hxb2_end']

                prefix = direction + '_'

                if 'X' in seq:
                    seqlen = len(seq)
                    seq = handle_x(seq)
                    if not seq or len(seq) < seqlen / 6:
                        new_row[prefix + 'error'] = errors.low_end_cov
                        continue

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
                                           sample_size=extended_length)
                    if finder2.is_full_length:
                        finder = finder2

                # Natalie's request
                # If a primer is not found at all, have a custom error for it, if there is something found but it did not pass secondary validation then make a different error for that
                if not finder.sample_primer:
                    new_row[prefix + 'error'] = errors.no_primer
                    continue
                elif not finder.is_valid:
                    new_row[prefix + 'error'] = errors.failed_validation
                    continue
                new_row[prefix +
                        'canonical_primer_seq'] = primers[direction]['seq']
                new_row[prefix + 'sample_primer_seq'] = finder.sample_primer
                new_row[prefix + 'sample_primer_start'] = finder.start
                new_row[prefix + 'sample_primer_end'] = finder.end
                new_row[prefix + 'sample_primer_size'] = len(
                    finder.sample_primer)
                new_row[prefix +
                        'hxb2_sample_primer_start'] = finder.hxb2_start
                new_row[prefix + 'hxb2_sample_primer_end'] = finder.hxb2_end

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


def load_csv(csv_filepath, name, seqtype, results=None):
    if results is None:
        results = {}
    df = pd.read_csv(csv_filepath)
    df['name'] = name
    df['seqtype'] = seqtype
    if name not in results:
        results[name] = {seqtype: df}
    else:
        results[name][seqtype] = df
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
        filtered = filtered.drop_duplicates(subset='sample', keep=False)
    # Remove any rows with references containing "reverse" or "unknown"
    filtered = filtered[(~filtered['reference'].str.contains('reverse'))
                        & (~filtered['reference'].str.contains('unknown'))]
    # duplicates = filtered.duplicated(subset='sample', keep=False)
    # duplicates = filtered[duplicates[duplicates].index]['sample'].unique()
    filtered = filtered[['name', 'sample', 'reference', 'sequence', 'seqtype']]
    return filtered


def run(contigs_csv,
        conseqs_csv,
        cascade_csv,
        name,
        outpath,
        disable_hivseqinr=False,
        nodups=True,
        split=1,
        sample_size=50,
        test=False):
    all_samples = utils.getSamplesFromCascade(cascade_csv)

    contigs_out = find_primers(contigs_csv,
                               outpath,
                               f'{name}_contigs',
                               seqtype='contigs',
                               sample_size=sample_size,
                               all_samples=all_samples,
                               test=test)
    conseqs_out = find_primers(conseqs_csv,
                               outpath,
                               f'{name}_conseqs',
                               seqtype='conseqs',
                               sample_size=sample_size,
                               all_samples=all_samples,
                               test=test)
    dfs = load_csv(contigs_out, name, 'contigs')
    dfs = load_csv(conseqs_out, name, 'conseqs', dfs)
    files = []
    for name in dfs:
        contigs_df = dfs[name]['contigs']
        conseqs_df = dfs[name]['conseqs']
        # Generate outcome summary
        OutcomeSummary(conseqs_df, contigs_df, outpath, test=test)
        # Generate the failure summary
        # utils.genFailureSummary(contigs_df, conseqs_df, outpath)
        filtered_contigs = filter_df(contigs_df, nodups)
        filtered_conseqs = filter_df(conseqs_df, nodups)
        joined = filtered_contigs.merge(filtered_conseqs,
                                        on='sample',
                                        suffixes=('_contig', '_conseq'),
                                        how='outer')
        joined.to_csv(outpath / 'joined.csv', index=False)
        joined['sequence'] = joined['sequence_conseq'].fillna(
            joined['sequence_contig'])
        joined['seqtype'] = joined['seqtype_conseq'].fillna(
            joined['seqtype_contig'])
        joined['name'] = joined['name_conseq'].fillna(joined['name_contig'])
        joined['seqlen'] = joined['sequence'].str.len()
        joined['reference'] = joined['reference_conseq'].fillna(
            joined['reference_contig'])
        joined = joined[[
            'name', 'sample', 'reference', 'seqtype', 'sequence', 'seqlen'
        ]].sort_values(by='sample')
        joined.to_csv(outpath / f'{name}_filtered.csv', index=False)
        nrows = len(joined)
        # How many rows per file
        nrows_per_file = math.ceil(nrows / split)
        for i in range(split):
            synthetic_primers_fasta = outpath / f'{name}_synthetic_primers_{i}.fasta'
            no_primers_fasta = outpath / f'{name}_no_primers_{i}.fasta'
            o = open(synthetic_primers_fasta, 'w')
            o2 = open(no_primers_fasta, 'w')
            start = i * nrows_per_file
            stop = start + nrows_per_file
            if stop > nrows:
                stop = nrows - 1
            for row in joined.iloc[start:stop].itertuples():
                # I don't remember why it was necessary to replace dashes with underscores but I think it was because HIVSEQINR doesn't like dashes in names
                # I've commented it out for now
                # header = f'>{row.name}_{row.sample}_{row.reference}_{row.seqtype}'.replace('-', '_')
                # The header delimiter, this must match the split in gene_splicer
                header = '>' + '::'.join(
                    (row.name, row.sample, row.reference, row.seqtype))
                o.write(
                    f'{header}\n{primers["fwd"]["nomix"] + row.sequence.replace("-", "") + primers["rev"]["nomix"]}\n'
                )
                o2.write(f'{header}\n{row.sequence.replace("-", "")}\n')
            o.close()
            o2.close()
            if not disable_hivseqinr:
                hivseqinr = Hivseqinr(outpath / f'hivseqinr_{i}',
                                      synthetic_primers_fasta)
            files.append(no_primers_fasta)
    return files


def main():
    args = parse_args()
    fasta_files = run(contigs_csv=args.contigs_csv,
                      conseqs_csv=args.conseqs_csv,
                      cascade_csv=args.cascade_csv,
                      name=args.name,
                      outpath=args.outpath.resolve(),
                      disable_hivseqinr=args.disable_hivseqinr,
                      nodups=args.nodups,
                      split=args.split,
                      sample_size=args.sample_size)
    return {'fasta_files': fasta_files, 'args': args}


if __name__ in ('__main__', '__live_coding__'):
    main()
