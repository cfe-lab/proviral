#! /usr/bin/env python

import logging
import argparse
from typing import Sequence
import runpy
import sys
from pathlib import Path
from importlib.metadata import version


PROGRAMS = ['sample', 'pipeline', 'study_summary', 'landscapes', 'hivseqinr']


def get_version() -> str:
    try:
        return str(version(__package__))
    except:
        return "development"


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run Proviral Pipeline.", add_help=False)
    parser.add_argument("--version", action="store_true", help="Print version and exit.")
    parser.add_argument('--help', action='store_true', help='Show this help message and exit.')
    parser.add_argument("program", nargs='?', choices=PROGRAMS, help="Program name.")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="Program arguments.")
    return parser


def execute_module_as_main(module_name: str, arguments: Sequence[str]) -> int:
    root_directory = str(Path(__file__).parent.parent)
    if root_directory not in sys.path:
        sys.path.append(root_directory)
    runpy.run_module(module_name, run_name='__main__')
    return 0


def execute_program(program: str, arguments: Sequence[str]) -> int:
    mod = f'cfeproviral.{program}'
    sys.argv = ['cfeproviral ' + program] + list(arguments)
    return execute_module_as_main(mod, arguments)


def main(argv: Sequence[str]) -> int:
    logging.basicConfig(level=logging.WARNING)
    parser = get_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(get_version())
        return 0
    elif args.help:
        parser.print_help()
        return 0
    elif args.program is None:
        parser.print_help()
        return 1

    program: str = args.program
    program_args: Sequence[str] = args.arguments
    return execute_program(program, program_args)


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__': entry()  # noqa
