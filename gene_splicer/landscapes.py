import csv
import logging
import os
import re
import typing
from typing import TextIO, Mapping, Dict, Set, List, Iterable, Tuple
import argparse
import sys

import json
import shutil
import subprocess as sp
import glob
from pathlib import Path
from csv import DictWriter, DictReader
from itertools import groupby
from operator import itemgetter

from gene_splicer.utils import (
    iterate_cfeintact_verdicts_1,
    iterate_hivseqinr_verdicts_1,
    LEFT_PRIMER_END, RIGHT_PRIMER_START,
)


logger = logging.getLogger(__name__)


def generate_proviral_landscape_csv_1_cont(blastn_reader: csv.DictReader,
                                           landscape_writer: csv.DictWriter,
                                           verdicts: Mapping[str, str],
                                           ) -> None:

    for row in blastn_reader:
        if row['qseqid'] in ['8E5LAV', 'HXB2']:
            # skip the positive control rows
            continue

        ref_start = int(row['sstart'])
        ref_end = int(row['send'])
        if ref_end <= LEFT_PRIMER_END or ref_start >= RIGHT_PRIMER_START:
            # skip unspecific matches of LTR at start and end
            continue

        qseqid = row['qseqid']
        try:
            [run_name, sample_name, _, _] = qseqid.split('::')
        except ValueError:
            [run_name, sample_name] = [None, qseqid]

        is_inverted = ''
        if ref_end < ref_start:
            # automatically recognize inverted regions
            new_end = ref_start
            ref_start = ref_end
            ref_end = new_end
            is_inverted = 'yes'

        if qseqid not in verdicts:
            logger.error("Could not generate landscapes for qseqid %r.", qseqid)
            continue

        verdict = verdicts[qseqid]
        is_defective = verdict != 'Intact'
        landscape_entry = {'ref_start': ref_start,
                           'ref_end': ref_end,
                           'samp_name': sample_name,
                           'run_name': run_name,
                           'is_inverted': is_inverted,
                           'is_defective': is_defective,
                           'defect': verdict,
                           }

        landscape_writer.writerow(landscape_entry)


def get_cfeintact_verdicts_1_map(details_dir: Path) -> Mapping[str, str]:
    ret: Dict[str, str] = {}

    for [qseqid, verdict] in iterate_cfeintact_verdicts_1(details_dir):
        ret[qseqid] = verdict

    return ret


def get_hivseqinr_verdicts_1_map(details_dir: Path) -> Mapping[str, str]:
    ret: Dict[str, str] = {}

    for [qseqid, verdict] in iterate_hivseqinr_verdicts_1(details_dir):
        ret[qseqid] = verdict

    return ret


def generate_proviral_landscape_csv_1(landscape_writer: csv.DictWriter,
                                      details_dir: Path,
                                      ) -> None:
    is_cfeintact = (details_dir / "holistic.csv").exists()
    if is_cfeintact:
        verdicts = get_cfeintact_verdicts_1_map(details_dir)
        blastn_path = details_dir / "blast.csv"
    else:
        verdicts = get_hivseqinr_verdicts_1_map(details_dir)
        blastn_path = details_dir / "Results_Intermediate" / "Output_Blastn_HXB2MEGA28_tabdelim.txt"

    with blastn_path.open() as blastn_file:
        if is_cfeintact:
            blastn_reader = DictReader(blastn_file)
        else:
            blastn_columns = ['qseqid',
                              'qlen',
                              'sseqid',
                              'sgi',
                              'slen',
                              'qstart',
                              'qend',
                              'sstart',
                              'send',
                              'evalue',
                              'bitscore',
                              'length',
                              'pident',
                              'nident',
                              'btop',
                              'stitle',
                              'sstrand']
            blastn_reader = DictReader(blastn_file, fieldnames=blastn_columns, delimiter='\t')

        return generate_proviral_landscape_csv_1_cont(
            blastn_reader,
            landscape_writer,
            verdicts,
        )


def generate_proviral_landscape_csv(outpath: Path, is_cfeintact: bool):
    proviral_landscape_csv = os.path.join(outpath, 'proviral_landscape.csv')

    if is_cfeintact:
        subpath = 'cfeintact*'
    else:
        subpath = 'hivseqinr*'

    landscape_columns = ['samp_name', 'run_name', 'ref_start', 'ref_end', 'defect', 'is_inverted', 'is_defective']
    with open(proviral_landscape_csv, 'w') as landscape_file:
        landscape_writer = csv.DictWriter(landscape_file, fieldnames=landscape_columns)
        landscape_writer.writeheader()

        for details_dir in outpath.glob(subpath):
            generate_proviral_landscape_csv_1(landscape_writer, details_dir)


class UserError(RuntimeError):
    def __init__(self, fmt: str, *fmt_args: object):
        self.fmt = fmt
        self.fmt_args = fmt_args
        self.code = 1


def dir_path(string: str) -> Path:
    if os.path.exists(string) and os.path.isdir(string):
        return Path(string)
    else:
        raise UserError("Path %r does not exist or is not a directory.", string)


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate Proviral Landscape CSV.")

    parser.add_argument("--details_dir", type=dir_path, required=True,
                        help="Directory containing details files for verdicts.")

    parser.add_argument("--output", type=argparse.FileType("wt"), required=True,
                        help="Output CSV file for proviral landscape.")

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true',
                                 help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true',
                                 help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true',
                                 help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true',
                                 help='Minimize output verbosity.')

    return parser


def main(argv: list) -> int:
    parser = cli_parser()
    args = parser.parse_args(argv)
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logger.debug("Start.")

    fieldnames = ['ref_start', 'ref_end', 'samp_name', 'run_name', 'is_inverted', 'is_defective', 'defect']

    landscape_writer = csv.DictWriter(args.output, fieldnames=fieldnames)
    landscape_writer.writeheader()
    generate_proviral_landscape_csv_1(landscape_writer, args.details_dir)

    logger.debug("Done.")
    return 0


if __name__ == '__main__':
    try:
        rc = main(sys.argv[1:])
    except BrokenPipeError:
        logger.debug("Broken pipe.")
        rc = 0
    except KeyboardInterrupt:
        logger.debug("Interrupted.")
        rc = 1
    except UserError as e:
        logger.fatal(e.fmt, *e.fmt_args)
        rc = e.code

    sys.exit(rc)
