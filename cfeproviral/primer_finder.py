import re
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
import os
from tarfile import TarFile
import cfeintact

import pandas as pd
from pathlib import Path
import logging
import math

import typing
from typing import Literal, Optional

from cfeproviral.primer_finder_errors import PrimerFinderErrors
from cfeproviral.primer_finder_class import PrimerFinder
from cfeproviral.outcome_summary import OutcomeSummary
from cfeproviral.hivseqinr import Hivseqinr
from cfeproviral.helpers.proviral_helper import ProviralHelper

import cfeproviral.utils as utils
mixture_dict = utils.mixture_dict
reverse_and_complement = utils.reverse_and_complement

logger = logging.getLogger(__name__)

ROOT = Path.cwd()
while ROOT.parent != ROOT:
    ROOT = ROOT.parent

HIVSEQINR_PATH = ROOT / "opt" / "hivseqinr"

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
        help='Length of sequence (probe) from each end of sample to search for '
             'primer',
        default=50)
    parser.add_argument('-o',
                        '--outpath',
                        help='The path to save the output',
                        default=Path(os.getcwd()),
                        type=Path)
    parser.add_argument('--hivseqinr',
                        action='store_true',
                        help="Launch the HIVSeqinR analysis.")
    parser.add_argument('--cfeintact',
                        action='store_true',
                        help="Launch the CFEIntact analysis.")
    parser.add_argument(
        '--nodups',
        action='store_false',
        help='Set this flag to disable the removal of duplicate samples that '
             'pass QC'
    )
    parser.add_argument(
        '--split',
        type=int,
        default=1,
        help='To avoid memory issues in hivseqinr, split the resulting '
             'qc-passed sequences into this number of fastas, each will be '
             'processed sequentially and then all will be merged into the '
             'final result'
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
        extended_length=250,
        force_all_proviral=False):
    proviral_helper = ProviralHelper(force_all_proviral=force_all_proviral)
    errors = PrimerFinderErrors()
    v3_reference = 'HIV1-CON-XX-Consensus-seed'
    make_path(outpath)
    columns = [
        'run_name', 'sample', 'reference', 'error', 'sequence', 'seqlen',
        'nmixtures', 'is_rev_comp'
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

    # all_samples is a dict whose keys are sample names and values are the remap
    # column (constructed from cascade.csv file)
    for sample in all_samples:
        # If no reads remapped, contig/conseq does not exist, write it as an error
        if all_samples[sample] == 0:
            # Do not analyze non-proviral samples
            if not proviral_helper.is_proviral(sample):
                logger.debug('Skipping sample "%s" because it is non-proviral' %
                             sample)
                continue
            new_row = dict(run_name=run_name,
                           sample=sample,
                           error=errors.no_sequence)
            writer.writerow(new_row)

    if 'sample' in reader.fieldnames:
        sample_groups = groupby(reader, itemgetter('sample'))
    else:
        sample_name, = all_samples.keys()  # Should be exactly one entry.
        sample_groups = [(sample_name, reader)]
    for sample_name, sample_rows in sample_groups:

        # Do not analyze non-proviral samples
        if not (force_all_proviral or proviral_helper.is_proviral(sample_name)):
            logger.debug('Skipping sample "%s" because it is non-proviral' %
                         sample_name)
            continue

        contig_num = 0
        non_hiv_rows = []

        for row in sample_rows:
            contig_num += 1
            seed_name = row.get('ref') or row['region'] or row.get('genotype')

            conseq_cutoff = row.get('consensus-percent-cutoff')
            contig_name = f'{contig_num}-{seed_name}'

            new_row = dict(run_name=run_name,
                           reference=contig_name,
                           is_rev_comp='N')
            new_row['sample'] = sample_name

            contig_seq: str = row.get('contig') or row['sequence']
            contig_seq = contig_seq.upper()

            if 'HIV' in seed_name:
                non_hiv_rows = None
            elif non_hiv_rows is not None:
                new_row['error'] = errors.non_hiv
                non_hiv_rows.append(new_row)
                continue
            else:
                continue

            # If percent consensus cutoff is not max, skip
            if conseq_cutoff and conseq_cutoff != 'MAX':
                new_row['error'] = errors.not_max
                writer.writerow(new_row)
                continue

            # Determine if sequence has internal Xs
            x_locations = [i for i, j in enumerate(contig_seq) if j == 'X']
            if any([(sample_size < i < len(contig_seq) - sample_size)
                    for i in x_locations]):
                new_row['error'] = errors.low_internal_cov
                writer.writerow(new_row)
                continue

            # Skip v3 sequences
            try:
                # If "region" is a column of row, then we are looking at a
                # conseq and not a contig. Only conseqs can have V3 sequences
                # so if we can't access this key we do nothing
                if row['region'] == v3_reference:
                    new_row['error'] = errors.non_proviral
                    writer.writerow(new_row)
                    continue
            except KeyError:
                pass

            new_row['seqlen'] = len(contig_seq)
            new_row['sequence'] = contig_seq

            # Determine if sequence has non-tcga characters
            found_non_tcga = re.findall(non_tcga, contig_seq)
            mixtures = len([x for x in found_non_tcga if x[0].upper() != 'X'])
            if mixtures > 1:
                new_row['error'] = errors.non_tcga
                new_row['nmixtures'] = mixtures
                writer.writerow(new_row)
                continue

            for key in columns:
                if key not in ('sample', 'contig', 'seqlen', 'error',
                               'sequence', 'run_name', 'reference',
                               'is_rev_comp'):
                    new_row[key] = None
            rev_row = dict(new_row)
            record_primers(contig_seq,
                           new_row,
                           errors,
                           sample_size,
                           extended_length)
            if new_row.get('error') is None and (
                    new_row.get('fwd_error') is not None or
                    new_row.get('rev_error') is not None):
                # Forward version failed, try reverse complement.
                utils.complement_dict['X'] = 'X'
                rev_seq = reverse_and_complement(contig_seq)
                del utils.complement_dict['X']
                record_primers(rev_seq,
                               rev_row,
                               errors,
                               sample_size,
                               extended_length)
                if (rev_row.get('error') is None and
                        rev_row.get('fwd_error') is None and
                        rev_row.get('rev_error') is None):
                    # Reverse complement worked, use it instead.
                    rev_row['is_rev_comp'] = 'Y'
                    new_row = rev_row
                    logger.warning(f'primers found in reverse complement of '
                                   f'{run_name}: {sample_name} - '
                                   f'{contig_name} - {seqtype}')

            writer.writerow(new_row)
        if non_hiv_rows is not None:
            if non_hiv_rows:
                writer.writerows(non_hiv_rows)
            elif all_samples[sample_name] != 0:
                new_row = dict(run_name=run_name,
                               sample=sample_name,
                               error=errors.no_sequence)
                writer.writerow(new_row)

    outfile.close()
    return outfilepath


def record_primers(contig_seq, new_row, errors, sample_size, extended_length):
    prime5_seq = contig_seq[:sample_size]
    extended_prime5_seq = contig_seq[:extended_length]
    prime3_seq = contig_seq[-sample_size:]
    extended_prime3_seq = contig_seq[-extended_length:]
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
        # If a primer is not found at all, have a custom error for it, if there
        # is something found but it did not pass secondary validation then make
        # a different error for that
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


def remove_primers(sample_size, row):
    # Strip the primers out, convert index values from floats.
    N = sample_size
    fwd_end = int(row.fwd_sample_primer_start + row.fwd_sample_primer_size)
    rev_start = len(row.sequence) - N + int(row.rev_sample_primer_start)
    if rev_start <= fwd_end:
        return row

    newseq = row.sequence[fwd_end:rev_start]
    row.sequence = newseq
    return row


def filter_df(sample_size, df, nodups=True):
    filtered = df[(df['error'].isna()
                   & df['fwd_error'].isna()
                   & df['rev_error'].isna())]
    filtered = filtered.apply(lambda x: remove_primers(sample_size, x), axis=1, result_type='broadcast')
    if nodups:
        filtered = filtered.drop_duplicates(subset='sample', keep=False)
    if not filtered.empty:
        # Remove any rows with references containing "reverse" or "unknown"
        filtered = filtered[(~filtered['reference'].str.contains('reverse'))
                            & (~filtered['reference'].str.contains('unknown'))]
    # duplicates = filtered.duplicated(subset='sample', keep=False)
    # duplicates = filtered[duplicates[duplicates].index]['sample'].unique()
    columns = ['reference', 'sequence', 'seqtype']
    if 'sample' in filtered.columns:
        columns.insert(0, 'sample')
    if 'name' in filtered.columns:
        columns.insert(0, 'name')
    filtered = filtered[columns]
    return filtered


def archive_hivseqinr_results(working_path: Path,
                              hivseqinr_results_tar: Path):
    final_results_path = working_path / 'Results_Final'
    with hivseqinr_results_tar.open("wb") as writer:
        archive = TarFile(fileobj=writer, mode='w')
        for result_path in final_results_path.iterdir():
            archive.add(result_path, result_path.name)


def archive_cfeintact_results(working_path: Path,
                              cfeintact_results_tar: Path):
    with cfeintact_results_tar.open("wb") as writer:
        archive = TarFile(fileobj=writer, mode='w')
        for result_path in working_path.iterdir():
            archive.add(result_path, result_path.name)


def run(contigs_csv,
        conseqs_csv,
        cascade_csv,
        name,
        outpath,
        nodups=True,
        split=1,
        sample_size=50,
        force_all_proviral=False,
        default_sample_name: Optional[str] = None,
        hivseqinr_results_tar: Optional[Path] = None,
        backend: utils.Backend = None,
        cfeintact_results_tar: Optional[Path] = None):
    all_samples = utils.get_samples_from_cascade(cascade_csv,
                                                 default_sample_name)

    contigs_out = find_primers(contigs_csv,
                               outpath,
                               f'{name}_contigs',
                               seqtype='contigs',
                               sample_size=sample_size,
                               all_samples=all_samples,
                               force_all_proviral=force_all_proviral)
    conseqs_out = find_primers(conseqs_csv,
                               outpath,
                               f'{name}_conseqs',
                               seqtype='conseqs',
                               sample_size=sample_size,
                               all_samples=all_samples,
                               force_all_proviral=force_all_proviral)
    dfs = load_csv(contigs_out, name, 'contigs')
    dfs = load_csv(conseqs_out, name, 'conseqs', dfs)
    files = []
    for name in dfs:
        contigs_df = dfs[name]['contigs']
        conseqs_df = dfs[name]['conseqs']
        # Generate outcome summary
        OutcomeSummary(conseqs_df, contigs_df, outpath, force_all_proviral)
        # Generate the failure summary
        # utils.genFailureSummary(contigs_df, conseqs_df, outpath)
        filtered_contigs = filter_df(sample_size, contigs_df, nodups)
        filtered_conseqs = filter_df(sample_size, conseqs_df, nodups)
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
        if joined.empty:
            joined['seqlen'] = 0
        else:
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
                # I don't remember why it was necessary to replace dashes with
                # underscores but I think it was because HIVSEQINR doesn't like dashes in names
                # I've commented it out for now
                # header = f'>{row.name}_{row.sample}_{row.reference}_{row.seqtype}'.replace('-', '_')
                # The header delimiter, this must match the split in gene_splicer
                header = f'>{row.name}::{row.sample}::{row.reference}::{row.seqtype}'
                o.write(
                    f'{header}\n{primers["fwd"]["nomix"] + row.sequence.replace("-", "") + primers["rev"]["nomix"]}\n'
                )
                o2.write(f'{header}\n{row.sequence.replace("-", "")}\n')
            o.close()
            o2.close()

            if backend == "HIVSeqinR":
                working_path = outpath / f'hivseqinr_{i}'
                hivseqinr_runner = Hivseqinr(HIVSEQINR_PATH,
                                             working_path,
                                             synthetic_primers_fasta)
                hivseqinr_runner.run()
                if hivseqinr_results_tar is not None:
                    archive_hivseqinr_results(working_path,
                                              hivseqinr_results_tar)

            elif backend == "CFEIntact":
                working_path = outpath / f'cfeintact_{i}'
                log_file_path = working_path / 'hiv-intact.log'
                os.makedirs(working_path, exist_ok=True)

                logger = cfeintact.logger
                file_handler = logging.FileHandler(log_file_path)
                logger.addHandler(file_handler)

                cfeintact.check(
                    output_dir=working_path,
                    input_file=str(no_primers_fasta),
                    subtype="B",
                    check_packaging_signal=True,
                    check_rre=True,
                    check_major_splice_donor_site=True,
                    check_hypermut=True,
                    check_long_deletion=True,
                    check_nonhiv=True,
                    check_scramble=True,
                    check_internal_inversion=True,
                    check_unknown_nucleotides=True,
                    check_small_orfs=True,
                    check_distance=False,
                    output_csv=True,
                )

                if cfeintact_results_tar is not None:
                    archive_cfeintact_results(working_path,
                                              cfeintact_results_tar)

            files.append(no_primers_fasta)
    return files


def main():
    args = parse_args()

    backend = None
    if args.cfeintact:
        backend = "CFEIntact"
    elif args.hivseqinr:
        backend = "HIVSeqinR"

    fasta_files = run(contigs_csv=args.contigs_csv,
                      conseqs_csv=args.conseqs_csv,
                      cascade_csv=args.cascade_csv,
                      name=args.name,
                      outpath=args.outpath.resolve(),
                      nodups=args.nodups,
                      split=args.split,
                      sample_size=args.sample_size,
                      backend=backend)

    return {'fasta_files': fasta_files, 'args': args}


if __name__ in ('__main__', '__live_coding__'):
    main()
