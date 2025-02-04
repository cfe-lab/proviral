#! /usr/bin/env python

import logging
import argparse
from typing import Sequence, Set, NoReturn
import runpy
import ast
import sys
import os
from io import StringIO
from pathlib import Path
import time
import traceback
import csv
import shlex


logger = logging.getLogger(__name__)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Start Proviral Pipeline daemon.")
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.', default=True)
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.')
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')
    return parser


ROOT = Path('.').absolute()
INPUT = ROOT / "input"


def process_directory(directory: Path) -> None:
    logger.info("Processing new directory: %r.", directory.name)
    os.chdir(directory)

    command_file = directory / "command"
    if command_file.exists():
        command_text = command_file.read_text()
    else:
        command_text = repr(["sample",

                             "--hivseqinr", "/opt/hivseqinr",

                             "sample_info.csv",
                             "contigs.csv",
                             "conseqs.csv",
                             "cascade.csv",

                             "outcome_summary.csv",
                             "conseqs_primers.csv",
                             "contigs_primers.csv",
                             "table_precursor.csv",
                             "proviral_landscape.csv",
                             "detailed_results.tar",
                             ])

    try:
        command: Sequence[str] = ast.literal_eval(command_text)
    except BaseException as ex:
        raise ValueError(f"Bad format of command, {command_text!r}.") from ex

    logger.info("Running %s.", ' '.join(map(shlex.quote, ['cfeproviral'] + list(command))))
    from gene_splicer.main import main
    try:
        ret = str(main(command))
    except SystemExit as ex:
        ret = str(ex)
    logger.info("Directory %r analyzed, exit code is %s.", directory.name, ret)


def try_process_directory(directory: Path) -> None:
    try:
        return process_directory(directory)
    except BaseException as ex:
        traceio = StringIO()
        traceback.print_exception(type(ex), ex, ex.__traceback__, file=traceio)
        traceio.seek(0)
        trace = traceio.read()
        logger.error("Failed to process directory %r:\n%s", directory.name, trace)


def serve() -> NoReturn:
    INPUT.mkdir(parents=True, exist_ok=True)
    last_directories: Set[str] = set()

    while True:
        current_directories = tuple(INPUT.iterdir())
        for directory in current_directories:
            if directory not in last_directories:
                try_process_directory(directory)

        last_directories = set(current_directories)
        time.sleep(0.1)


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logger.info("Starting the proviral daemon...")
    logger.info("Files copied to proviral:input will be processed by the Pipeline.")
    serve()


def entry() -> None:
    try:
        ret = main(sys.argv[1:])
    except KeyboardInterrupt:
        ret = 1
    sys.exit(ret)


if __name__ == '__main__': entry()  # noqa
