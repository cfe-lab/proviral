import os
import sys
import yaml
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
inputs = cwd / 'inputs'
outputs = cwd / 'outputs'
sys.path.append(str(cwd.parent.parent))

import utils

utils.write_fasta({'MOD_HXB2': utils.mod_hxb2}, inputs / 'mod_hxb2.fasta')


def test_intact():
    # For this test I will simply use the full HXB2 reference as my sequence, this should guarantee it's intact-ness
    alignment = utils.load_samfile(inputs / 'alignment.sam')
    coords = utils.splice_genes(utils.mod_hxb2, utils.mod_hxb2, alignment,
                                utils.mod_annot)
    genes = utils.coords_to_genes(coords, utils.mod_hxb2)
    valid_genes = utils.load_yaml(inputs / 'valid_genes.yaml')
    for gene, seq in genes.items():
        assert seq == valid_genes[gene]