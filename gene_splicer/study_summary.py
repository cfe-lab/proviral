import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictWriter, DictReader

import typing
from io import TextIOBase
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from subprocess import run, CalledProcessError, PIPE, STDOUT


def parse_args():
    parser = ArgumentParser(
        description='Run HIVSeqinR on many run directories, then summarize.',
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog="This will create a subfolder for each run in the same folder "
               "as the study summary CSV file. If a run's subfolder already "
               "exists and has an output_summary.csv file, then the run won't "
               "be reprocessed.")
    parser.add_argument('study_summary_csv',
                        default='study_summary.csv',
                        nargs='?',
                        type=FileType('w'),
                        help='Study summary, plus a path for run output folders.')
    parser.add_argument('--samples_csv',
                        type=FileType(),
                        help='CSV file that has a run column in YYMMDD_MACHINE '
                             'format, as well as sample and pid (participant '
                             'id) columns.')
    parser.add_argument('--runs_root',
                        default='.',
                        type=Path,
                        help='Root folder to find runs listed in samples_csv.')
    parser.add_argument('--runs',
                        nargs='*',
                        help='Choose which runs to summarize. Defaults to all '
                             'the runs in samples_csv. Any samples not found '
                             'in samples.csv will guess the participant id '
                             'from the first part of the sample name.')
    return parser.parse_args()


class StudySummary:
    def __init__(self):
        self.run_paths: typing.Tuple[Path, ...] = ()

        # {run_path: {col_name: count}}
        self.run_counts: typing.Dict[str, typing.Dict[str, int]] = \
            defaultdict(Counter)

        # {participant_id: {col_name: count}}
        self.participant_counts: typing.Dict[str, typing.Dict[str, int]] = \
            defaultdict(Counter)

        # {(run_name, sample): participant_id}
        self.participant_ids: typing.Dict[typing.Tuple[str, str], str] = {}

        # [(run_name, sample)]
        self.unmapped_samples: typing.List[typing.Tuple[str, str]] = []

    def load_runs(self, requested_runs: typing.Iterable[str]):
        self.run_paths = tuple(Path(run_path) for run_path in requested_runs)

    def load_samples(self, samples_csv: TextIOBase, runs_root: Path):
        old_run_paths = set(self.run_paths)
        run_selections = {}  # {run_name: is_selected}
        new_run_paths = set()
        for row in DictReader(samples_csv):
            run_name = row['run']
            is_selected = run_selections.get(run_name)
            if is_selected is None:
                is_selected = False
                for run_path in runs_root.glob(run_name + '*'):
                    if old_run_paths and run_path not in old_run_paths:
                        continue
                    is_selected = True
                    new_run_paths.add(run_path)
                run_selections[run_name] = is_selected
            if is_selected:
                self.participant_ids[(run_name, row['sample'])] = row['pid']
        if not old_run_paths:
            self.run_paths = tuple(sorted(new_run_paths))

    def load_outcome(self, outcome_summary_csv: typing.TextIO):
        for row in DictReader(outcome_summary_csv):
            run_name = row['run']
            run_counts = self.run_counts[run_name]
            sample_name = row['sample']
            participant_id_guess = sample_name.split('-')[0]
            sample_key = (run_name, sample_name)
            participant_id = self.participant_ids.get(sample_key)
            if participant_id is None:
                participant_id = participant_id_guess
                self.unmapped_samples.append(sample_key)
            participant_counts = self.participant_counts[participant_id]

            for counts in (run_counts, participant_counts):
                counts['samples'] += 1
                passed_field = row['passed']
                if passed_field == 'True':
                    counts['passed'] += 1
                else:
                    assert passed_field == 'False', passed_field
                error_type = row['error']
                if error_type:
                    counts['errors'] += 1
                    column_name = error_type.replace(' ', '_')
                    counts[column_name] += 1

    def write(self, study_summary_csv: typing.TextIO):
        count_columns = ['samples',
                         'passed',
                         'errors',
                         'no_sequence',
                         'non_hiv',
                         'primer_error',
                         'low_internal_cov',
                         'multiple_contigs']
        writer = DictWriter(study_summary_csv,
                            ['type',  # run, participant, or total
                             'name'] + count_columns,
                            lineterminator=os.linesep)
        writer.writeheader()
        total_counts = Counter()
        for run_name, run_counts in sorted(self.run_counts.items()):
            total_counts += run_counts
            for column_name in count_columns:
                run_counts[column_name] += 0  # Force zeroes in CSV.
            row = dict(type='run', name=run_name, **run_counts)
            writer.writerow(row)

        for participant_id, participant_counts in sorted(
                self.participant_counts.items()):
            for column_name in count_columns:
                participant_counts[column_name] += 0  # Force zeroes in CSV.
            row = dict(type='participant',
                       name=participant_id,
                       **participant_counts)
            writer.writerow(row)

        for column_name in count_columns:
            total_counts[column_name] += 0  # Force zeroes in CSV.
        total_row = dict(type='total', name='total', **total_counts)
        writer.writerow(total_row)

    def write_warnings(self, report_file: typing.TextIO, limit: int = None):
        if not self.unmapped_samples:
            return
        print('WARNING, some samples did not map to participant ids:',
              file=report_file)
        sample_count = 0
        skipped_sample_count = 0
        skipped_run_count = 0
        for run_name, items in groupby(self.unmapped_samples, itemgetter(0)):
            if limit is not None and sample_count > limit:
                skipped_run_count += 1
                skipped_sample_count += sum(1 for _ in items)
                continue
            print(run_name + ':', file=report_file)
            for _, sample in items:
                if limit is not None:
                    sample_count += 1
                    if sample_count > limit:
                        skipped_count = 1 + sum(1 for _ in items)
                        print(f'  [{skipped_count} more samples...]',
                              file=report_file)
                        break
                print(' ', sample, file=report_file)
        if skipped_run_count:
            print(f'[{skipped_sample_count} more samples in '
                  f'{skipped_run_count} more runs...]',
                  file=report_file)


def run_gene_splicer(run_path: Path, outcome_folder: Path):
    version_results_path = run_path / 'Results' / 'version_7.14'
    assert version_results_path.exists(), version_results_path
    denovo_path = version_results_path / 'denovo'
    contigs_path = str(denovo_path / 'contigs.csv')
    conseq_path = str(denovo_path / 'conseq.csv')
    cascade_path = str(denovo_path / 'cascade.csv')
    outcome_folder.mkdir(exist_ok=True)
    full_run_name = outcome_folder.name
    short_run_name = '_'.join(full_run_name.split('_')[:2])
    python_path = sys.executable
    launch_path = Path(__file__).parent.parent
    log_path = outcome_folder / 'gene_splicer.log'
    pipeline_args = [python_path,
                     '-m', 'gene_splicer.pipeline',
                     '--outpath', str(outcome_folder),
                     contigs_path,
                     conseq_path,
                     cascade_path,
                     short_run_name]
    try:
        with log_path.open('w') as log_file:
            run(pipeline_args,
                stdout=log_file,
                stderr=STDOUT,
                check=True,
                encoding=sys.getdefaultencoding(),
                cwd=launch_path)
    except CalledProcessError:
        print(log_path.read_text(), file=sys.stderr)
        raise
    print('.', end='')


def main():
    args = parse_args()
    summary = StudySummary()
    if args.runs:
        summary.load_runs(args.runs)
    if args.samples_csv:
        summary.load_samples(args.samples_csv, args.runs_root)

    outcome_root = Path(args.study_summary_csv.name).parent
    launch_count = 0
    for run_path in summary.run_paths:
        run_folder_name = run_path.name
        outcome_path: Path = (outcome_root / run_folder_name /
                              'outcome_summary.csv')
        if not outcome_path.exists():
            run_gene_splicer(run_path, outcome_path.parent)
            assert outcome_path.exists(), outcome_path
            launch_count += 1
        with outcome_path.open() as outcome_summary_csv:
            try:
                summary.load_outcome(outcome_summary_csv)
            except Exception as ex:
                raise RuntimeError(
                    f'Failed to load outcome from {outcome_path}.') from ex

    if launch_count:
        print()  # New line after ... progress display.

    with args.study_summary_csv:
        summary.write(args.study_summary_csv)
    summary.write_warnings(sys.stdout, limit=100)


if __name__ == '__main__':
    main()
