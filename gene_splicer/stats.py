import argparse
import os
import sys
import csv
import pandas
import yaml
from pathlib import Path
from typing import ContextManager, List
from contextlib import contextmanager
from gene_splicer.primer_finder_errors import PrimerFinderErrors
from gene_splicer.helpers.proviral_helper import ProviralHelper

errors = PrimerFinderErrors()


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute statistics for proviral runs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--contigs_analysis_csvs',
                        required=True,
                        nargs='+',
                        help='A csv file produced by the proviral pipeline')
    parser.add_argument('--filtered_csvs',
                        required=True,
                        nargs='+',
                        help='A csv file produced by the proviral pipeline')
    parser.add_argument(
        '--sample_mapping',
        type=Path,
        help='A csv file with 2 columns: sample and participant_id')
    parser.add_argument('--outpath',
                        type=Path,
                        default=Path(os.getcwd()).resolve(),
                        help='Path to output files')
    parser.add_argument(
        '--force_all_proviral',
        action='store_true',
        help=
        'FOR TESTING PURPOSES. Forces all samples to be considered proviral.')
    args = parser.parse_args()
    helper = ProviralHelper(force_all_proviral=args.force_all_proviral)
    if not args.sample_mapping:
        args.sample_mapping = helper.get_sample_pid_mapping()
    return args


class OutputFile:
    def __init__(self, outpath: Path, name: str, ftype: str = 'csv') -> None:
        self.outpath = outpath
        self.name = name
        self.ftype = ftype
        self.path = self.get_path()

    def get_path(self):
        return self.outpath / (self.name + '.' + self.ftype)


class CsvFile(OutputFile):
    def __init__(self,
                 outpath: Path,
                 name: str,
                 fieldnames: list = None,
                 ftype: str = 'csv') -> None:
        super().__init__(outpath, name, ftype)
        self.fieldnames = fieldnames

    def init_write(self):
        self.csvfile = open(self.path, 'w', newline='')
        self.writer = csv.DictWriter(self.csvfile, fieldnames=self.fieldnames)
        if self.fieldnames:
            self.writer.writeheader()

    def write_rows(self, data):
        with open(self, 'w') as o:
            for item in data:
                self.writer.writerow(item)

    def close(self):
        self.csvfile.close()

    @contextmanager
    def open(self):
        self.init_write()
        try:
            yield self.writer
        finally:
            self.csvfile.close()


def get_unique_samples_from_contigs_analysis(contigs_analysis_csv_path):
    csv = pandas.read_csv(contigs_analysis_csv_path)
    # return set(csv[csv['sample'].str.contains('HIV')]['sample'].unique())
    return set(csv['sample'].unique())


def main():
    args = parse_args()
    run(args=args)


def get_unique_samples(contigs_files):
    unique_samples = set()
    for csv in contigs_files:
        unique_samples = unique_samples.union(
            get_unique_samples_from_contigs_analysis(csv))
    unique_samples_data = []
    fieldnames = ['sample_name', 'is_proviral']
    for unique_sample in unique_samples:
        is_proviral = False
        if 'HIV' in unique_sample:
            is_proviral = True
        row = {'sample_name': unique_sample, 'is_proviral': is_proviral}
        unique_samples_data.append(row)
    unique_samples_data.sort(key=lambda x: x['is_proviral'])
    unique_samples_outfile = CsvFile('.',
                                     'unique_samples_data',
                                     fieldnames=fieldnames)
    unique_samples_outfile.write_rows(unique_samples_data)


def init_run(row):
    if row['error'] == errors.no_sequence:
        return errors.no_sequence
    elif 'unknown' in row['reference'] or row['error'] in (errors.is_v3,
                                                           errors.non_hiv):
        return errors.non_hiv
    elif row['fwd_error'] or row['rev_error']:
        return errors.hiv_but_failed
    else:
        return True


def get_data(contigs_file, filtered_file, data=None, force_all_proviral=False):
    proviral_helper = ProviralHelper(force_all_proviral=force_all_proviral)
    if not data:
        data = {}
    df = pandas.read_csv(contigs_file)
    run = None
    for index, row in df.iterrows():
        if not proviral_helper.is_proviral(row['sample']):
            continue
        sample = row['sample']
        run = row['run_name']
        if sample not in data:
            data[sample] = {run: init_run(row)}
        elif run not in data[sample]:
            data[sample][run] = init_run(row)

    df = pandas.read_csv(filtered_file)
    for index, row in df.iterrows():
        sample = row['sample']
        if not proviral_helper.is_proviral(sample):
            continue
        # Check if sample has already passed, meaning two contigs passed
        if data[sample][run] is True:
            print(f'WARNING: Run {run} already passed sample {sample}!')
        data[sample][run] = True
    return data


def compute_failure_rate(subset, total):
    try:
        return round(100 * (subset / total))
    except ZeroDivisionError:
        return 0


def compute_per_participant(data, mapping_filepath):
    mapping = {}
    with open(mapping_filepath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mapping[row['sample']] = row['participant_id']
    result = {}
    for sample, runs in data.items():
        pid = mapping[sample]
        if pid not in result:
            result[pid] = {
                errors.no_sequence: 0,
                errors.non_hiv: 0,
                errors.hiv_but_failed: 0,
                'passed': 0,
                'total': 0
            }
        for run, verdict in runs.items():
            if verdict == errors.no_sequence:
                result[pid][errors.no_sequence] += 1
            elif verdict == errors.non_hiv:
                result[pid][errors.non_hiv] += 1
            elif verdict == errors.hiv_but_failed:
                result[pid][errors.hiv_but_failed] += 1
            elif verdict == 'passed':
                result[pid]['passed'] += 1
            result[pid]['total'] += 1
    for pid, qcs in result.items():
        failed = qcs['total'] - qcs['passed']
        percent_failure = compute_failure_rate(failed, qcs['total'])
        percent_failure_no_sequence = compute_failure_rate(
            qcs[errors.no_sequence], failed)
        percent_failure_non_hiv = compute_failure_rate(qcs[errors.non_hiv],
                                                       failed)
        percent_failure_hiv = compute_failure_rate(qcs[errors.hiv_but_failed],
                                                   failed)
        result[pid] = {
            'total': qcs['total'],
            'passed': qcs['passed'],
            'failed': failed,
            'percent_failure': percent_failure,
            'failed_no_sequence': qcs[errors.no_sequence],
            'percent_failed_no_sequence': percent_failure_no_sequence,
            'failed_non_hiv': qcs[errors.non_hiv],
            'percent_failed_non_hiv': percent_failure_non_hiv,
            'failed_hiv': qcs[errors.hiv_but_failed],
            'percent_failed_hiv': percent_failure_hiv
        }
    return result


def compute_stats(data):
    per_run = {}
    total = len(data)
    total_passed = 0
    for sample, runs in data.items():
        for run, verdict in runs.items():
            if run not in per_run:
                per_run[run] = {
                    'total': 0,
                    'passed': 0,
                    errors.no_sequence: 0,
                    errors.non_hiv: 0,
                    errors.hiv_but_failed: 0,
                    'percent_failure': {}
                }
            per_run[run]['total'] += 1
            if verdict is True:
                per_run[run]['passed'] += 1
                total_passed += 1
            elif verdict == errors.no_sequence:
                per_run[run][errors.no_sequence] += 1
            elif verdict == errors.non_hiv:
                per_run[run][errors.non_hiv] += 1
            elif verdict == errors.hiv_but_failed:
                per_run[run][errors.hiv_but_failed] += 1
    for run in per_run:
        per_run[run]['failed'] = per_run[run]['total'] - per_run[run]['passed']
        per_run[run]['total_percent_failure'] = compute_failure_rate(
            per_run[run]['failed'], per_run[run]['total'])
        # Get the percentage of failures that are due to no sequence
        per_run[run]['percent_failure'][
            errors.no_sequence] = compute_failure_rate(
                per_run[run][errors.no_sequence], per_run[run]['failed'])
        # Get the percentage of failures that are due to non hiv
        per_run[run]['percent_failure'][errors.non_hiv] = compute_failure_rate(
            per_run[run][errors.non_hiv], per_run[run]['failed'])
        # Get the percentage of failures that are due to some other failure
        per_run[run]['percent_failure'][
            errors.hiv_but_failed] = compute_failure_rate(
                per_run[run][errors.hiv_but_failed], per_run[run]['failed'])
    overall_percent_failure = compute_failure_rate(total - total_passed, total)
    return overall_percent_failure, per_run


def get_all_data(contigs_files_paths,
                 filtered_files_paths,
                 force_all_proviral=False):
    assert len(contigs_files_paths) == len(filtered_files_paths)
    i = 0
    all_data = {}
    while i < len(contigs_files_paths):
        contig_path = contigs_files_paths[i]
        filtered_path = filtered_files_paths[i]
        all_data = get_data(contig_path,
                            filtered_path,
                            all_data,
                            force_all_proviral=force_all_proviral)
        i += 1
    return all_data


def write_data(data, outpath):
    output = CsvFile(outpath,
                     'cumulative_data',
                     fieldnames=['run', 'sample', 'passed'])
    with output.open() as o:
        for sample, runs in data.items():
            for run, passed in runs.items():
                o.writerow({'run': run, 'sample': sample, 'passed': passed})


def write_stats(per_run, overall, outpath):
    output = CsvFile(outpath,
                     'per_run',
                     fieldnames=[
                         'run', 'total', 'passed', 'failed', 'percent_failure',
                         'failed_no_sequence', 'percent_failed_no_sequence',
                         'failed_non_hiv', 'percent_failed_non_hiv',
                         'failed_hiv', 'percent_failed_hiv'
                     ])
    with output.open() as o:
        for run, value in per_run.items():
            o.writerow({
                'run':
                run,
                'failed':
                value['failed'],
                'passed':
                value['passed'],
                'total':
                value['total'],
                'percent_failure':
                value['total_percent_failure'],
                'failed_no_sequence':
                value[errors.no_sequence],
                'percent_failed_no_sequence':
                value['percent_failure'][errors.no_sequence],
                'failed_non_hiv':
                value[errors.non_hiv],
                'percent_failed_non_hiv':
                value['percent_failure'][errors.non_hiv],
                'failed_hiv':
                value[errors.hiv_but_failed],
                'percent_failed_hiv':
                value['percent_failure'][errors.hiv_but_failed]
            })
    print(f'Overall failure rate: {overall}')


def write_per_participant(data, outpath):
    output = CsvFile(outpath,
                     'per_participant',
                     fieldnames=[
                         'pid', 'total', 'passed', 'failed', 'percent_failure',
                         'failed_no_sequence', 'percent_failed_no_sequence',
                         'failed_non_hiv', 'percent_failed_non_hiv',
                         'failed_hiv', 'percent_failed_hiv'
                     ])
    with output.open() as o:
        for pid, mystats in data.items():
            o.writerow({
                'pid':
                pid,
                'failed':
                mystats['failed'],
                'passed':
                mystats['passed'],
                'total':
                mystats['total'],
                'percent_failure':
                mystats['percent_failure'],
                'failed_no_sequence':
                mystats['failed_no_sequence'],
                'percent_failed_no_sequence':
                mystats['percent_failed_no_sequence'],
                'failed_non_hiv':
                mystats['failed_non_hiv'],
                'percent_failed_non_hiv':
                mystats['percent_failed_non_hiv'],
                'failed_hiv':
                mystats['failed_hiv'],
                'percent_failed_hiv':
                mystats['percent_failed_hiv']
            })


def run(args):
    # unique_samples = get_unique_samples(args.contigs_analysis_csvs)
    all_data = get_all_data(args.contigs_analysis_csvs,
                            args.filtered_csvs,
                            force_all_proviral=args.force_all_proviral)
    write_data(all_data, args.outpath)
    overall_percent_failure, per_run = compute_stats(all_data)
    write_stats(per_run, overall_percent_failure, args.outpath)
    per_participant = compute_per_participant(all_data, args.sample_mapping)
    write_per_participant(per_participant, args.outpath)


if __name__ == '__main__':
    main()
