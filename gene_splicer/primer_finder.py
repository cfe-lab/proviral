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
import logger
from primer_finder_class import PrimerFinder

import utils
mixture_dict = utils.mixture_dict
reverse_and_complement = utils.reverse_and_complement

# from micall.core.project_config import ProjectConfig
# from micall.utils.translation import mixture_dict, reverse_and_complement
from probe_finder import ProbeFinder

logger = logging.getLogger('gene_splicer')


class Hivseqinr:
    def __init__(self, outpath, fasta):
        self.url = 'https://github.com/guineverelee/HIVSeqinR/raw/master/HIVSeqinR_ver2.7.1.zip'
        self.outpath = outpath
        self.fasta = fasta
        self.download()
        self.make_blast_dir()
        self.fix_wds()
        self.copy_fasta()
        self.job = self.run()
        self.finalize()

    def to_rpath(self, path):
        return path.replace('\\', '\\\\')

    def copy_fasta(self):
        raw_fastas_path = self.outpath / 'RAW_FASTA'
        try:
            os.makedirs(raw_fastas_path)
        except FileExistsError:
            pass
        logger.debug('Attempting to copy "%s" to "%s"' %
                     (self.fasta, raw_fastas_path))
        shutil.copy(self.fasta, raw_fastas_path)
        return True

    def make_blast_dir(self):
        dbdir = self.outpath / 'hxb2_blast_db'
        try:
            os.mkdir(dbdir)
        except FileExistsError:
            pass
        hxb2_path = dbdir / 'HXB2.fasta'
        shutil.copyfile(self.outpath / 'R_HXB2.fasta', hxb2_path)
        cmd = [
            'makeblastdb', '-in', hxb2_path, '-parse_seqids', '-dbtype', 'nucl'
        ]
        self.dbdir = dbdir
        job = subprocess.run(cmd)
        return job

    def download(self):
        response = requests.get(self.url)
        file_like_object = io.BytesIO(response.content)
        zipfile_obj = zipfile.ZipFile(file_like_object)
        if not os.path.isdir(self.outpath):
            os.makedirs(self.outpath)
        zipfile_obj.extractall(self.outpath)
        return True

    def fix_wds(self):
        rscript_path = os.path.join(
            self.outpath, 'R_HIVSeqinR_Combined_ver09_ScrambleFix.R')
        new_rscript_path = os.path.join(self.outpath, 'modified.R')
        self.new_rscript_path = new_rscript_path
        with open(self.new_rscript_path, 'w') as outfile:
            with open(rscript_path, 'r') as infile:
                for line in infile:
                    if 'MyWD <- getwd()' in line:
                        line = line.replace('MyWD <- getwd()',
                                            f'MyWD = "{self.outpath}"\n')
                    elif line.startswith('MyBlastnDir <-'):
                        line = f'MyBlastnDir = "{str(self.dbdir) + os.path.sep*2}"\n'
                    outfile.write(line)

    def run(self):
        cwd = os.getcwd()
        os.chdir(self.outpath)
        cmd = ['Rscript', './modified.R']
        job = subprocess.run(cmd)
        os.chdir(cwd)
        return job

    def finalize(self):
        path = os.path.join(self.outpath, 'COMPLETE')
        Path(path).touch()


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


def find_primers(csv_filepath,
                 outpath,
                 run_name,
                 sample_size=50,
                 extended_length=200):
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
    skipped = {}
    total = 0
    viable = 0
    unique_samples = 0
    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):
        # Do not analyze non-proviral samples
        if not utils.is_proviral(sample_name):
            logger.debug('Skipping sample "%s" because it is non-proviral' %
                         sample_name)
            continue
        contig_num = 0
        unique_samples += 1
        for row in sample_rows:
            total += 1
            seed_name = row.get('genotype') or row.get('ref') or row['region']
            conseq_cutoff = row.get('consensus-percent-cutoff')
            contig_num += 1
            contig_name = f'{contig_num}_{seed_name}'
            uname = f'{sample_name}_{contig_name}_{contig_num}'

            # interesting_sample = 'HIV3428G2-P19-HIV_S56'
            # pause = False
            # if ((sample_name == interesting_sample)
            #  and (seed_name == '1-HIV1-B-FR-K03455-seed')
            #  and (run_name.endswith('conseqs'))):
            #     pause = True

            new_row = dict(run_name=run_name,
                           sample=sample_name,
                           reference=contig_name)
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

            try:
                # If "region" is a column of row, then we are looking at a conseq and not a contig. Only conseqs can have V3 sequences so if we can't access this key we do nothing
                if row['region'] == v3_reference:
                    skipped[uname] = 'is V3 sequence'
                    new_row['error'] = skipped[uname]
                    writer.writerow(new_row)
                    continue
            except KeyError:
                pass

            # interest = 'HIV3428F1-L22-HIV_S6'
            # if (sample_name == interest
            #     and contig_name == '1_1-HIV1-B-FR-K03455-seed'
            # ):
            #     import pdb; pdb.set_trace()

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
            extended_prime5_seq = contig_seq[:extended_length]
            prime3_seq = contig_seq[-sample_size:]
            extended_prime3_seq = contig_seq[-extended_length:]
            for key in columns:
                if key not in [
                        'sample', 'contig', 'seqlen', 'error', 'sequence',
                        'run_name', 'reference'
                ]:
                    new_row[key] = None

            # interest = 'T0288D90-M7-HIV_S95'
            # if (sample_name == interest
            #         and contig_name == '1_1-HIV1-B-FR-K03455-seed'):
            #     import pdb
            #     pdb.set_trace()

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
                        skipped[uname] = 'too many X in sequence'
                        new_row[prefix + 'error'] = skipped[uname]
                        continue

                # interest = 'T0288D90-M7-HIV_S95'
                # if (sample_name == interest
                #         and contig_name == '8_2-unknown-partial'):
                #     import pdb
                #     pdb.set_trace()

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

                # This can happen if there are too many combinations of the sequence to unpack (too many mixtures or dashes etc)
                if not finder.is_valid:
                    skipped[uname] = 'primer was not found'
                    new_row[prefix + 'error'] = skipped[uname]
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
            viable += 1
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


def run(contigs_csv, conseqs_csv, name, outpath, disable_hivseqinr, nodups,
        split, sample_size):
    contigs_out = find_primers(contigs_csv, outpath, f'{name}_contigs',
                               sample_size)
    conseqs_out = find_primers(conseqs_csv, outpath, f'{name}_conseqs',
                               sample_size)
    dfs = load_csv(contigs_out, name, 'contigs')
    dfs = load_csv(conseqs_out, name, 'conseqs', dfs)
    files = []
    for name in dfs:
        contigs_df = dfs[name]['contigs']
        conseqs_df = dfs[name]['conseqs']
        # Generate the failure summary
        utils.genFailureSummary(contigs_df, conseqs_df, outpath)
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
                      name=args.name,
                      outpath=args.outpath.resolve(),
                      disable_hivseqinr=args.disable_hivseqinr,
                      nodups=args.nodups,
                      split=args.split,
                      sample_size=args.sample_size)
    return {'fasta_files': fasta_files, 'args': args}


if __name__ in ('__main__', '__live_coding__'):
    main()
