import argparse
import os
import csv
from pathlib import Path
from cfeproviral.primer_finder_errors import PrimerFinderErrors
from cfeproviral.helpers.proviral_helper import ProviralHelper


# Base this file on the outcome summary. For all samples total, for all samples in each run, and for all samples belonging to every participant compute the following columns:
# total, passed, each of the 5 possible failure types in the outcome summary "error" column
def parse_args():
    parser = argparse.ArgumentParser(
        description='Compute statistics for proviral runs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outcome_summary_csv',
                        required=True,
                        nargs='+',
                        help='A csv file produced by the proviral pipeline')
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
    return args


def main():
    args = parse_args()
    run(args=args)


def run(args):
    errors = PrimerFinderErrors()
    pvhelper = ProviralHelper(sample_map_path_override=args.sample_mapping,
                              force_all_proviral=args.force_all_proviral)


if __name__ == '__main__':
    main()
