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

import gene_splicer.utils as utils

utils.write_fasta({'MOD_HXB2': utils.mod_hxb2}, cwd / 'mod_hxb2.fasta')
valid_genes = utils.load_yaml(cwd / 'valid' / 'valid_genes.yaml')


def test_large_deletion1():
    name = 'large_deletion1'
    # This should just be GAG now
    query = utils.mod_hxb2[123:1626]
    utils.write_fasta({'GAG': query}, cwd / f'{name}.fasta')
    # I also verified this sequence by using LANL sequence locator
    gag = valid_genes['GAG']
    # Just sanity check that my query sequence is GAG
    assert gag == query
    # Align the query to the target
    # I have to skip this section because minimap2 does not run on windows
    # Otherwise I could call it as a subprocess
    # For now I will manually align the query to the target
    alignment = utils.load_samfile(cwd / f'{name}.sam')
    # Get the genes
    results = utils.splice_genes(query, utils.mod_hxb2, alignment,
                                 utils.mod_annot)
    genes = utils.coords_to_genes(results, query)
    for gene in valid_genes:
        try:
            assert genes[gene] == valid_genes[gene]
        except AssertionError as e:
            print(f'Gene {gene} does not match')
            raise AssertionError(e)
