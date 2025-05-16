import os
from pathlib import Path

import pandas as pd

import cfeproviral.primer_finder as primer_finder

cwd = Path(os.path.realpath(__file__)).parent
inputs = cwd.parent / 'data' / 'example3' / 'inputs'
outputs = cwd.parent / 'data' / 'example3' / 'outputs'


def test_filter_df():
    sample_size = 50
    sample = pd.DataFrame(
        data={
            'error': [''],
            'fwd_error': [''],
            'rev_error': [''],
            'reference': [''],
            'sequence': [''],
            'seqtype': ['']
        })
    primer_finder.filter_df(sample_size, sample)

    # This breaks because of a string filter being applied to a numeric column
    # sample = pd.DataFrame(
    #     data={
    #         'error': [''],
    #         'fwd_error': [''],
    #         'rev_error': [''],
    #         'reference': [5],
    #         'sequence': [''],
    #         'seqtype': ['']
    #     })
    # filtered = primer_finder.filter_df(sample)

    # this used to return an error when NaN values were in the reference column
    # and pandas was trying to do a str operation!
    contigs_csv_path = inputs / 'contigs.csv'
    conseqs_csv_path = inputs / 'conseqs.csv'
    outpath = outputs
    all_samples = {'my_sample': 0}
    run_name = 'my_run'
    contigs_out = primer_finder.find_primers(open(contigs_csv_path),
                                             outpath,
                                             run_name,
                                             all_samples,
                                             'contigs',
                                             force_all_proviral=True)
    conseqs_out = primer_finder.find_primers(open(conseqs_csv_path),
                                             outpath,
                                             run_name,
                                             all_samples,
                                             'conseqs',
                                             force_all_proviral=True)
    dfs = primer_finder.load_csv(contigs_out, run_name, 'contigs')
    dfs = primer_finder.load_csv(conseqs_out, run_name, 'conseqs', dfs)
    # If you remove the below 2 lines, the NaN (which are actually floats) will
    # break the filter_df function since it uses a pandas .str filter internally
    contigs_df = dfs[run_name]['contigs'].fillna('')
    conseqs_df = dfs[run_name]['conseqs'].fillna('')
    primer_finder.filter_df(sample_size, contigs_df, nodups=True)
    primer_finder.filter_df(sample_size, conseqs_df, nodups=True)
