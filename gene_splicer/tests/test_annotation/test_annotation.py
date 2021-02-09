import os
import sys
import subprocess
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

import gene_splicer.utils as utils


def test_annotation():
    actual_annot = utils.mod_annot
    expected_annot = {x['gene']: [int(x['start']), int(x['stop'])] for x in utils.read_csv(cwd / 'valid' / 'valid_mod_annot.csv')}
    utils.csv_to_bed(cwd / 'valid' / 'valid_mod_annot.csv', 'MOD_HXB2')
    assert actual_annot == expected_annot