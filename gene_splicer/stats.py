import argparse
import os
import sys
import csv
import pandas
from pathlib import Path
from typing import ContextManager, List
from contextlib import contextmanager


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
        required=True,
        help='A csv file with 2 columns: sample and participant_id')
    parser.add_argument('--outpath',
                        type=Path,
                        default=Path(os.getcwd()).resolve(),
                        help='Path to output files')
    return parser.parse_args()


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


def isProviral(sample_name):
    if (('GAGGAG' in sample_name) or ('VIR' in sample_name) or
        ('NEF-HIV') in sample_name) or ('V3LOOP' in sample_name) or (
            'HLA' in sample_name) or ('HCV' in sample_name):
        return False
    return True


def get_data(contigs_file, filtered_file, data=None):
    if not data:
        data = {}
    df = pandas.read_csv(contigs_file)
    run = None
    for index, row in df.iterrows():
        if not isProviral(row['sample']):
            continue
        sample = row['sample']
        run = row['run_name']
        if sample not in data:
            data[sample] = {run: False}
        elif run not in data[sample]:
            data[sample][run] = False

    df = pandas.read_csv(filtered_file)
    for index, row in df.iterrows():
        sample = row['sample']
        if not isProviral(sample):
            continue
        # Check if sample has already passed, meaning two contigs passed
        if data[sample][run] is True:
            print(f'WARNING: Run {run} already passed sample {sample}!')
        data[sample][run] = True
    return data


def compute_failure_rate(total, passed):
    return round(100 * ((total - passed) / total))


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
            result[pid] = []
        for run, passed in runs.items():
            result[pid].append(passed)
    for pid, qcs in result.items():
        total = len(qcs)
        passed = 0
        for qc in qcs:
            if qc:
                passed += 1
        failed = total - passed
        percent_failure = compute_failure_rate(total, passed)
        result[pid] = {
            'total': total,
            'passed': passed,
            'failed': failed,
            'percent_failure': percent_failure
        }
    return result


def compute_stats(data):
    per_run = {}
    total = len(data)
    total_passed = 0
    for sample, runs in data.items():
        for run, passed in runs.items():
            if run not in per_run:
                per_run[run] = {'total': 0, 'passed': 0}
            per_run[run]['total'] += 1
            if passed:
                per_run[run]['passed'] += 1
                total_passed += 1
    for run in per_run:
        per_run[run]['failed'] = per_run[run]['total'] - per_run[run]['passed']
        per_run[run]['percent_failure'] = compute_failure_rate(
            per_run[run]['total'], per_run[run]['passed'])
    overall_percent_failure = compute_failure_rate(total, total_passed)
    return overall_percent_failure, per_run


def get_all_data(contigs_files_paths, filtered_files_paths):
    assert len(contigs_files_paths) == len(filtered_files_paths)
    i = 0
    all_data = {}
    while i < len(contigs_files_paths):
        contig_path = contigs_files_paths[i]
        filtered_path = filtered_files_paths[i]
        all_data = get_data(contig_path, filtered_path, all_data)
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
    output = CsvFile(
        outpath,
        'per_run',
        fieldnames=['run', 'total', 'passed', 'failed', 'percent_failure'])
    with output.open() as o:
        for run, value in per_run.items():
            o.writerow({
                'run': run,
                'failed': value['failed'],
                'passed': value['passed'],
                'total': value['total'],
                'percent_failure': value['percent_failure']
            })
    print(f'Overall failure rate: {overall}')


def write_per_participant(data, outpath):
    output = CsvFile(
        outpath,
        'per_participant',
        fieldnames=['pid', 'total', 'passed', 'failed', 'percent_failure'])
    with output.open() as o:
        for pid, mystats in data.items():
            o.writerow({
                'pid': pid,
                'failed': mystats['failed'],
                'passed': mystats['passed'],
                'total': mystats['total'],
                'percent_failure': mystats['percent_failure']
            })


def run(args):
    # unique_samples = get_unique_samples(args.contigs_analysis_csvs)
    all_data = get_all_data(args.contigs_analysis_csvs, args.filtered_csvs)
    write_data(all_data, args.outpath)
    overall_percent_failure, per_run = compute_stats(all_data)
    write_stats(per_run, overall_percent_failure, args.outpath)
    per_participant = compute_per_participant(all_data, args.sample_mapping)
    write_per_participant(per_participant, args.outpath)


if __name__ == '__main__':
    main()