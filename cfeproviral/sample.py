import logging
import typing
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader
from pathlib import Path
from importlib.metadata import version
from typing import Sequence

import cfeproviral.gene_splicer as gene_splicer
import cfeproviral.primer_finder as primer_finder
import cfeproviral.utils as utils
import cfeproviral.landscapes as landscapes


def output_file(argument: str) -> Path:
    path = Path(argument)
    if (not path.exists()) or path.is_file():
        return path
    else:
        raise ValueError(f"Path {argument!r} is not a regular file.")


def get_argument_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description='Search sequences from a single sample for primers, and '
                    'identify errors.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    # inputs
    parser.add_argument('sample_info_csv',
                        help='sample name and other details',
                        type=FileType())
    parser.add_argument('contigs_csv',
                        help='CSV file with contig sequences',
                        type=FileType())
    parser.add_argument('conseqs_csv',
                        help='CSV file with conseq sequences',
                        type=FileType())
    parser.add_argument('cascade_csv',
                        help='CSV file from MiCall output',
                        type=FileType())
    # outputs
    parser.add_argument('outcome_summary_csv',
                        help='Summary result for the whole sample',
                        type=output_file)
    parser.add_argument('conseqs_primers_csv',
                        help='Analysis of primers in consensus sequences',
                        type=output_file)
    parser.add_argument('contigs_primers_csv',
                        help='Analysis of primers in assembled contig sequences',
                        type=output_file)
    parser.add_argument('table_precursor_csv',
                        help='Sequence data ready to upload',
                        type=output_file)
    parser.add_argument('proviral_landscape_csv',
                        help='Data for proviral landscape plot',
                        type=output_file)
    parser.add_argument('detailed_results_tar',
                        help="Archive file with HIVSeqinR's final results "
                             "folder, or CFEIntact's results.",
                        type=output_file)
    parser.add_argument(
        '-p',
        '--sample_size',
        type=int,
        help='Length of sequence (probe) from each end of sample to search for primer',
        default=50)
    parser.add_argument('--hivseqinr',
                        action='store_true',
                        help="Launch the HIVSeqinR analysis.")
    parser.add_argument('--cfeintact',
                        action='store_true',
                        help="Launch the CFEIntact analysis.")
    parser.add_argument(
        '--nodups',
        action='store_false',
        help='Set this flag to disable the removal of duplicate samples that pass QC'
    )
    parser.add_argument(
        '--split',
        type=int,
        default=1,
        help='To avoid memory issues in hivseqinr, split the resulting '
             'qc-passed sequences into this number of fastas, each will be '
             'processed sequentially and then all will be merged into the '
             'final result. Obsolete for CFEIntact.')
    return parser


def copy_output(source: Path, target: Path):
    source.rename(target)


def main(argv: Sequence[str]) -> int:
    logging.basicConfig(level=logging.WARNING)
    parser = get_argument_parser()
    args = parser.parse_args(argv)
    outpath = args.outcome_summary_csv.parent / 'scratch'
    outpath = outpath.resolve()
    outpath.mkdir(exist_ok=True, parents=True)

    # Read sample_info to get run_name, sample, and build the mapping
    info_reader = DictReader(args.sample_info_csv)
    rows = list(info_reader)

    if rows:
        first_row = rows[0]
        run_name = first_row.get('run_name', 'kive_run')
        default_sample_name = first_row.get('sample')  # None if no 'sample' column

        # Build sample_info_map
        if info_reader.fieldnames and 'sample' in info_reader.fieldnames:
            # Multi-sample format: each row maps to a specific sample
            sample_info_map = {row['sample']: row for row in rows}
            if not default_sample_name:
                default_sample_name = rows[0]['sample']  # Use first sample as default
        else:
            # No 'sample' column: this is run-level metadata that applies to all samples
            # Leave sample_info_map empty - we'll pass run_level_info separately
            sample_info_map = {}
            if not default_sample_name:
                default_sample_name = 'Sample_S1'  # Fallback default
    else:
        run_name = 'kive_run'
        default_sample_name = 'Sample_S1'
        sample_info_map = {}
        first_row = {}

    hivseqinr_results_tar = None
    cfeintact_results_tar = None
    backend: utils.Backend = None

    if args.cfeintact:
        cfeintact_results_tar = args.detailed_results_tar
        backend = "CFEIntact"
    elif args.hivseqinr:
        hivseqinr_results_tar = args.detailed_results_tar
        backend = "HIVSeqinR"

    fasta_files = primer_finder.run(contigs_csv=args.contigs_csv,
                                    conseqs_csv=args.conseqs_csv,
                                    cascade_csv=args.cascade_csv,
                                    hivseqinr_results_tar=hivseqinr_results_tar,
                                    name=run_name,
                                    outpath=outpath,
                                    nodups=args.nodups,
                                    split=args.split,
                                    sample_size=args.sample_size,
                                    force_all_proviral=True,
                                    default_sample_name=default_sample_name,
                                    backend=backend,
                                    cfeintact_results_tar=cfeintact_results_tar,
                                    sample_info_map=sample_info_map,
                                    run_level_info=first_row if not sample_info_map else None)
    for file in fasta_files:
        gene_splicer.run(file, outdir=outpath)
    utils.generate_table_precursor(name=run_name, outpath=outpath)
    landscapes.generate_proviral_landscape_csv(outpath, backend)
    copy_output(outpath / 'outcome_summary.csv', args.outcome_summary_csv)
    copy_output(outpath / (run_name + '_conseqs_primer_analysis.csv'),
                args.conseqs_primers_csv)
    copy_output(outpath / (run_name + '_contigs_primer_analysis.csv'),
                args.contigs_primers_csv)
    copy_output(outpath / 'table_precursor.csv', args.table_precursor_csv)
    copy_output(outpath / 'proviral_landscape.csv', args.proviral_landscape_csv)
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__': entry()  # noqa
