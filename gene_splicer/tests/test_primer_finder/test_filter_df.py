import os
import sys
import yaml
from pathlib import Path
import pandas as pd
import io

cwd = Path(os.path.realpath(__file__)).parent
inputs = cwd / 'inputs'
outputs = cwd / 'outputs'

import gene_splicer.primer_finder as primer_finder


def test_filter_df():
    sample = pd.DataFrame(
        data={
            'error': [''],
            'fwd_error': [''],
            'rev_error': [''],
            'reference': [''],
            'sequence': [''],
            'seqtype': ['']
        })
    filtered = primer_finder.filter_df(sample)

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

    # this used to return an error when NaN values were in the reference column and pandas was trying to do a str operation!
    contigs_csv_path = inputs / 'contigs.csv'
    conseqs_csv_path = inputs / 'conseqs.csv'
    outpath = outputs
    contigs_out = primer_finder.find_primers(open(contigs_csv_path), outpath,
                                             'contigs')
    conseqs_out = primer_finder.find_primers(open(conseqs_csv_path), outpath,
                                             'conseqs')
    dfs = primer_finder.load_csv(contigs_out, 'contigs')
    dfs = primer_finder.load_csv(conseqs_out, 'conseqs', dfs)
    files = []
    contigs_df = dfs['contigs'].fillna('')
    conseqs_df = dfs['conseqs'].fillna('')
    filtered_contigs = primer_finder.filter_df(contigs_df, nodups=True)
    filtered_conseqs = primer_finder.filter_df(conseqs_df, nodups=True)