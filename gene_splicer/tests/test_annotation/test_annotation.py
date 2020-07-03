import os
import sys
import subprocess
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

# This is a hack to get the tests running
# Unfortunately in order to use one of the python dependencies I need actual
# software downloaded and available via the command-line
# Therefore a singularity image would be the way to go, but it does not make
# for rapid testing since I have to have the singularity image working
# in order to run any tests. For now I will test this way to test the actual
# code and then use special sets of inputs/outputs to sort of do an
# integration on the singularity image
sys.path.append(str(cwd.parent.parent))

import utils


def test_annotation():
    actual_annot = utils.mod_annot
    expected_annot = {x['gene']: [int(x['start']), int(x['stop'])] for x in utils.read_csv(cwd / 'valid' / 'valid_mod_annot.csv')}
    utils.csv_to_bed(cwd / 'valid' / 'valid_mod_annot.csv', 'MOD_HXB2')
    assert actual_annot == expected_annot