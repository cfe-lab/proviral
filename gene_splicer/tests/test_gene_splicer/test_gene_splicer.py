import os
import sys
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


valid_genes = utils.load_yaml(cwd / 'valid' / 'valid_genes.yaml')


def test_large_deletion1():
    # This should just be GAG now
    large_del = utils.mod_hxb2[123:1626]
    # I also verified this sequence by using LANL sequence locator
    gag = valid_genes['GAG']
    assert gag == large_del