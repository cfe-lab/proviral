import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import Counter, defaultdict
from csv import DictWriter, DictReader

import typing
from io import TextIOBase
from pathlib import Path


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
        self.run_paths = tuple(sorted(new_run_paths))

    def load_outcome(self, outcome_summary_csv: typing.TextIO):
        for row in DictReader(outcome_summary_csv):
            run_name = row['run']
            run_counts = self.run_counts[run_name]
            sample_name = row['sample']
            participant_id_guess = sample_name.split('-')[0]
            participant_id = self.participant_ids.get((run_name, sample_name),
                                                      participant_id_guess)
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
                         'hiv_but_failed',
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


def main():
    args = parse_args()
    summary = StudySummary()
    if args.runs:
        summary.load_runs(args.runs)
    if args.samples_csv:
        summary.load_samples(args.samples_csv, args.runs_root)

    outcome_root = Path(args.study_summary_csv.name).parent
    for run_path in summary.run_paths:
        run_folder_name = run_path.name
        outcome_path: Path = (outcome_root / run_folder_name /
                              'outcome_summary.csv')
        if not outcome_path.exists():
            print(run_folder_name, 'missing.')
            continue
        with outcome_path.open() as outcome_summary_csv:
            summary.load_outcome(outcome_summary_csv)

    with args.study_summary_csv:
        summary.write(args.study_summary_csv)


if __name__ == '__main__':
    main()
