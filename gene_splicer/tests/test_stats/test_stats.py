import pytest
import os
import csv
from pathlib import Path

import gene_splicer.stats as stats

cwd = Path(os.path.realpath(__file__)).parent


class Args:
    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            setattr(self, k, v)


def test_stats():
    args = Args(contigs_analysis_csvs=[
        cwd / 'input' / 'contigs_primer_analysis_mod.csv'
    ],
                filtered_csvs=[cwd / 'input' / 'filtered_mod.csv'],
                outpath=(cwd / 'output'),
                sample_mapping=(cwd / 'input' / 'sample_mapping_mod.csv'),
                force_all_proviral=True)
    stats.run(args)
