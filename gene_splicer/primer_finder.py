import re
import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter
import logging
import os
import pandas as pd
import subprocess
import requests
import io
import zipfile
import shutil
from pathlib import Path
import sys

import Levenshtein
from gotoh import align_it

import utils
mixture_dict = utils.mixture_dict
reverse_and_complement = utils.reverse_and_complement

# from micall.core.project_config import ProjectConfig
# from micall.utils.translation import mixture_dict, reverse_and_complement
from probe_finder import ProbeFinder

logger = logging.getLogger('micall')

class Hivseqinr:


    def __init__(self, outpath, fasta):
        self.url = 'https://github.com/guineverelee/HIVSeqinR/raw/master/HIVSeqinR_ver2.7.1.zip'
        self.outpath = outpath
        self.fasta = fasta
        self.download()
        self.make_blast_dir()
        self.fix_wds()
        self.copy_fasta()
        self.run()
        self.finalize()


    def to_rpath(self, path):
        return path.replace('\\', '\\\\')


    def copy_fasta(self):
        raw_fastas_path = os.path.join(self.outpath, 'RAW_FASTA')
        os.mkdir(raw_fastas_path)
        shutil.copy(self.fasta, raw_fastas_path)
        return True


    def make_blast_dir(self):
        dbdir = os.path.join(self.outpath, 'hxb2_blast_db')
        os.mkdir(dbdir)
        hxb2_path = os.path.join(dbdir, 'HXB2.fasta')
        shutil.copyfile(
            os.path.join(self.outpath, 'R_HXB2.fasta'),
            hxb2_path
        )
        cmd = [
            'makeblastdb',
            '-in',
            hxb2_path,
            '-parse_seqids',
            '-dbtype',
            'nucl'
        ]
        self.dbdir = dbdir
        job = subprocess.run(cmd)
        return job


    def download(self):
        response = requests.get(self.url)
        file_like_object = io.BytesIO(response.content)
        zipfile_obj = zipfile.ZipFile(file_like_object)
        if not os.path.isdir(self.outpath):
            os.mkdir(self.outpath)
        zipfile_obj.extractall(self.outpath)
        return True


    def fix_wds(self):
        rscript_path = os.path.join(self.outpath, 'R_HIVSeqinR_Combined_ver09_ScrambleFix.R')
        new_rscript_path = os.path.join(self.outpath, 'modified.R')
        self.new_rscript_path = new_rscript_path
        with open(self.new_rscript_path, 'w') as outfile:
            with open(rscript_path, 'r') as infile:
                for line in infile:
                    if 'MyWD <- getwd()' in line:
                        line = line.replace('MyWD <- getwd()', f'MyWD = "{self.to_rpath(self.outpath)}"\n')
                    elif line.startswith('MyBlastnDir <-'):
                        line = f'MyBlastnDir = "{self.to_rpath(self.dbdir) + os.path.sep*2}"\n'
                    outfile.write(line)


    def run(self):
        cwd = os.getcwd()
        os.chdir(self.outpath)
        cmd = [
            'Rscript',
            './modified.R'
        ]
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
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences',
                        type=FileType())
    parser.add_argument('conseqs_csv',
                        help='CSV file with conseq sequences',
                        type=FileType())
    parser.add_argument('name',
                        help='A name for the analysis')
    parser.add_argument('-o',
                        '--outpath',
                        help='The path to save the output',
                        default=Path(os.getcwd()),
                        type=Path)
    parser.add_argument('--disable_hivseqinr',
                        action='store_true',
                        help='Disable running hivseqinr')
    return parser.parse_args()


def find_primers(csv_filepath, outpath, run_name):
    print(run_name)
    columns = ['run_name', 'sample', 'reference', 'error', 'sequence', 'seqlen', 'nmixtures']
    for target_name in primers:
        for column_type in [
            'probe_hxb2_start',
            'full_real_primer_seq',
            'in_probe_start',
            'in_probe_size',
            'in_hxb2_start',
            'in_hxb2_size',
            'is_reversed',
            'seq',
            'actual_primer_seq',
            'overhang',
            'finder_dist',
            'dist',
            'error',
            'warning',
        ]:
            columns.append(target_name + '_' + column_type)
    non_tcga = re.compile(r'[^TCGA-]+')
    outfilepath = os.path.join(outpath, f'{run_name}_primer_analysis.csv')
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
    unique_samples = 0
    for sample_name, sample_rows in groupby(reader, itemgetter('sample')):
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

            new_row = dict(run_name=run_name, sample=sample_name, reference=contig_name)
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

            # interest = 'HIV3428F1-L22-HIV_S6'
            # if (sample_name == interest
            #     and contig_name == '1_1-HIV1-B-FR-K03455-seed'
            # ):
            #     import pdb; pdb.set_trace()

            # Determine if sequence has internal Xs
            x_locations = [i for i,j in enumerate(contig_seq) if j=='X']
            probelen = 30
            if any([(probelen < i < len(contig_seq) - (probelen)) for i in x_locations]):
                skipped[uname] = 'contig sequence contained internal X'
                new_row['error'] = skipped[uname]
                writer.writerow(new_row)
                continue
            found_non_tcga = re.findall(non_tcga, contig_seq)
            mixtures = len([x for x in found_non_tcga if x[0].upper() != 'X'])
            if (
                mixtures > 1
            ):
                skipped[uname] = 'contig sequence contained non-TCGA/gap'
                new_row['error'] = skipped[uname]
                new_row['nmixtures'] = mixtures
                writer.writerow(new_row)
                continue
            prime5_seq = contig_seq[:probelen]
            prime3_seq = contig_seq[-probelen:]
            gap_open_penalty = 15
            gap_extend_penalty = 3
            use_terminal_gap_penalty = 1
            for key in columns:
                if key not in ['sample', 'contig', 'seqlen', 'error', 'sequence', 'run_name', 'reference']:
                    new_row[key] = None
            for end, seq in [(5, prime5_seq), (3, prime3_seq)]:
                if end == 5:
                    name = 'fwd'
                    hxb2_target_start = primers[name]['hxb2_start']
                    hxb2_target_end = primers[name]['hxb2_end'] + 100
                    hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
                else:
                    name = 'rev'
                    hxb2_target_start = primers[name]['hxb2_start'] - 100
                    hxb2_target_end = primers[name]['hxb2_end']
                    hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
                primer = None

                if 'X' in seq:
                    seqlen = len(seq)
                    oldseq = seq
                    seq = handle_x(seq)
                    if not seq or len(seq) < seqlen / 6:
                        skipped[uname] = 'too many X in sequence'
                        new_row[prefix + 'error'] = skipped[uname]
                        continue

                # if 'X' in seq:
                #     import pdb; pdb.set_trace()

                finder = ProbeFinder(hxb2_target_seq, seq)
                finder.start += hxb2_target_start
                prefix = name + '_'
                new_row[prefix + 'probe_hxb2_start'] = finder.start
                new_row[prefix + 'full_real_primer_seq'] = primers[name]['seq']

                # If the segment overlaps the primer
                if primers[name]['hxb2_start'] - probelen <= finder.start <= primers[name]['hxb2_end']:
                    primer = validate_primer(finder, seq, primers[name], hxb2_target_seq)
                    if primer['error']:
                        skipped[uname] = primer['error']
                        new_row[prefix + 'error'] = skipped[uname]
                # Otherwise
                else:
                    # If contig ends before hxb2 primer start
                    if primers[name]['hxb2_start'] - probelen > finder.start:
                        skipped[uname] = f'{end} contig probe ends before hxb2 primer start'
                    # If contig starts after hxb2 primer start
                    elif finder.start > primers[name]['hxb2_end']:
                        skipped[uname] = f'{end} contig probe starts after hxb2 primer end'
                    new_row[prefix + 'error'] = skipped[uname]
                    primer = validate_primer(finder, seq, primers[name], hxb2_target_seq)
                    if primer['error']:
                        skipped[uname] = primer['error']
                        new_row[prefix + 'error'] += f' primer_error: {skipped[uname]}'

                new_row[prefix + 'finder_dist'] = primer['finder_dist']
                if primer['finder_dist'] > (probelen / 10):
                    new_row[prefix + 'warning'] = f'Finder distance > {probelen / 10}'

                new_row[prefix + 'in_probe_start'] = primer['seq_start']
                try:
                    new_row[prefix + 'in_probe_size'] = primer['seq_end'] - primer['seq_start']
                except TypeError:
                    new_row[prefix + 'in_probe_size'] = None
                new_row[prefix + 'in_hxb2_start'] = primer['hxb2_start']
                try:
                    new_row[prefix + 'in_hxb2_size'] = primer['hxb2_end'] - primer['hxb2_start']
                except TypeError:
                    new_row[prefix + 'in_hxb2_size'] = None
                new_row[prefix + 'is_reversed'] = ('Y'
                                                    if finder.is_reversed
                                                    else 'N')
                new_row[prefix + 'seq'] = primer['target_seq']
                new_row[prefix + 'overhang'] = primer['overhang']
                new_row[prefix + 'dist'] = primer['dist']
                new_row[prefix + 'actual_primer_seq'] = primer['real_primer']
            writer.writerow(new_row)
            viable += 1
    outfile.close()
    return outfilepath


def handle_x(sequence):
    length = len(sequence)
    midpoint = length / 2
    x_positions = [i for i,j in enumerate(sequence) if j == 'X']
    try:
        rightmost_x = max([x for x in x_positions if x <= midpoint])
    except ValueError:
        rightmost_x = None
    try:
        leftmost_x = min([x for x in x_positions if x > midpoint])
    except ValueError:
        leftmost_x = None
    if rightmost_x is not None and leftmost_x is not None:
        sequence = sequence[rightmost_x+1:leftmost_x]
    elif rightmost_x is not None:
        sequence = sequence[rightmost_x+1:]
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
    primer_in_finder_hxb2_end_coord = min(
        finder.start + finder_seqsize,
        target['hxb2_end']
    )
    # Get the coordinates of the primer relative to the finder sequence
    primer_in_finder_start_coord = primer_in_finder_hxb2_start_coord - finder.start
    primer_in_finder_end_coord = primer_in_finder_hxb2_end_coord - finder.start
    if (
        (target['hxb2_start'] == 9603)
        and (primer_in_finder_end_coord >= matched_finder_size - overhang)
    ):
        primer_in_finder_start_coord -= overhang
        primer_in_finder_end_coord -= overhang

    # Get the primer sequence of the finder sequence
    primer_in_finder = finder_seq[primer_in_finder_start_coord:primer_in_finder_end_coord]

    # Get the sequence of the true primer that overlaps the finder sequence
    primer_start_coord = max(
        0,
        primer_in_finder_hxb2_start_coord - target['hxb2_start']
    )
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
        result['error'] = 'Real primer did not match length of primer in finder'
        return result
    if not real_primer:
        result['error'] = 'real primer not found at expected coordinates'
        return result
    elif not primer_in_finder:
        result['error'] = 'primer in contig sequence not found at expected coordinates'
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


def load_csv(csv_filepath, name, seqtype, results=None):
    if results is None:
        results = {}
    df = pd.read_csv(csv_filepath)
    df['name'] = name
    df['seqtype'] = seqtype
    if name not in results:
        results[name] = {
            seqtype: df
        }
    else:
        results[name][seqtype] = df
    return results


def add_primers(row):
    # Strip the primers out
    newseq = row.sequence[int(row.fwd_in_probe_size):-int(row.rev_in_probe_size)]
    # Add the primers in
#     newseq = primers['fwd']['seq'] + newseq + primers['rev']['seq']
    row.sequence = newseq
    return row


def remove_primers(row):
    # Strip the primers out
    newseq = row.sequence[int(row.fwd_in_probe_size):-int(row.rev_in_probe_size)]
    row.sequence = newseq
    return row


def filter_df(df):
    filtered = df[(
        df['error'].isna()
        & df['fwd_error'].isna()
        & df['rev_error'].isna()
    )]
    filtered = filtered.apply(remove_primers, axis=1)
    # filtered = filtered.drop_duplicates(subset='sample', keep=False)
    # duplicates = filtered.duplicated(subset='sample', keep=False)
    # duplicates = filtered[duplicates[duplicates].index]['sample'].unique()
    filtered = filtered[['name', 'sample', 'reference', 'sequence', 'seqtype']]
    return filtered


def output_filtered_data(contigs_csv, conseqs_csv, name, outpath, disable_hivseqinr):
    contigs_out = find_primers(contigs_csv, outpath, f'{name}_contigs')
    conseqs_out = find_primers(conseqs_csv, outpath, f'{name}_conseqs')
    dfs = load_csv(contigs_out, name, 'contigs')
    dfs = load_csv(conseqs_out, name, 'conseqs', dfs)
    for name in dfs:
        contigs_df = dfs[name]['contigs']
        conseqs_df = dfs[name]['conseqs']
        filtered_contigs = filter_df(contigs_df)
        filtered_conseqs = filter_df(conseqs_df)
        joined = filtered_contigs.merge(
            filtered_conseqs,
            on='sample',
            suffixes=('_contig', '_conseq'),
            how='outer'
        )
        joined['sequence'] = joined['sequence_conseq'].fillna(joined['sequence_contig'])
        joined['seqtype'] = joined['seqtype_conseq'].fillna(joined['seqtype_contig'])
        joined['name'] = joined['name_conseq'].fillna(joined['name_contig'])
        joined['seqlen'] = joined['sequence'].str.len()
        joined['reference'] = joined['reference_conseq'].fillna(joined['reference_contig'])
        joined = joined[[
            'name',
            'sample',
            'reference',
            'seqtype',
            'sequence',
            'seqlen'
        ]].sort_values(
            by='sample'
        )
        joined.to_csv(
            os.path.join(outpath, f'{name}_filtered.csv'),
            index=False
        )
        fasta_outpath = os.path.join(outpath, f'{name}_synthetic_primers.fasta')
        fasta_outpath2 = Path(outpath) / f'{name}_no_primers.fasta'
        o = open(fasta_outpath, 'w')
        o2 = open(fasta_outpath2, 'w')
        for row in joined.itertuples():
            # I don't remember why it was necessary to replace dashes with underscores but I think it was because HIVSEQINR doesn't like dashes in names
            header = f'>{row.name}_{row.sample}_{row.reference}_{row.seqtype}'.replace('-', '_')
            header = '>' + '::'.join((
                row.name,
                row.sample,
                row.reference,
                row.seqtype
            ))
            o.write(f'{header}\n{primers["fwd"]["nomix"] + row.sequence.replace("-", "") + primers["rev"]["nomix"]}\n')
            o2.write(f'{header}\n{row.sequence.replace("-", "")}\n')
        if not disable_hivseqinr:
            hivseqinr = Hivseqinr(os.path.join(outpath, 'hivseqinr'), fasta_outpath)
        o.close()
        o2.close()


def main():
    args = parse_args()
    output_filtered_data(args.contigs_csv, args.conseqs_csv, args.name, args.outpath, args.disable_hivseqinr)


if __name__ in ('__main__', '__live_coding__'):
    main()