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
                        type=FileType('w'))
    parser.add_argument('conseqs_primers_csv',
                        help='Analysis of primers in consensus sequences',
                        type=FileType('w'))
    parser.add_argument('contigs_primers_csv',
                        help='Analysis of primers in assembled contig sequences',
                        type=FileType('w'))
    parser.add_argument('table_precursor_csv',
                        help='Sequence data ready to upload',
                        type=FileType('w'))
    parser.add_argument('proviral_landscape_csv',
                        help='Data for proviral landscape plot',
                        type=FileType('w'))
    parser.add_argument('detailed_results_tar',
                        help="Archive file with HIVSeqinR's final results "
                             "folder, or CFEIntact's results.",
                        type=FileType('wb'))
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


def copy_output(source: Path, target: typing.IO):
    target.close()
    source.rename(target.name)


def main(argv: Sequence[str]) -> int:
    logging.basicConfig(level=logging.WARNING)
    parser = get_argument_parser()
    args = parser.parse_args(argv)
    outpath = Path(args.outcome_summary_csv.name).parent / 'scratch'
    outpath = outpath.resolve()
    outpath.mkdir(exist_ok=True, parents=True)
    with args.sample_info_csv:
        info_reader = DictReader(args.sample_info_csv)
        for row in info_reader:
            sample_info: dict = row
            break
        else:
            sample_info = {'run_name': 'test_run', 'sample': 'Sample_S1'}
    run_name = sample_info.get('run_name', 'kive_run')

    hivseqinr_results_tar = None
    cfeintact_results_tar = None
    backend = None

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
                                    default_sample_name=sample_info['sample'],
                                    backend=backend,
                                    cfeintact_results_tar=cfeintact_results_tar)
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
