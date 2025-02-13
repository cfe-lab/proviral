import os
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

import cfeproviral.utils as utils

utils.write_fasta({'MOD_HXB2': utils.mod_hxb2}, cwd / 'mod_hxb2.fasta')
valid_genes = utils.load_yaml(cwd / 'valid' / 'valid_genes.yaml')


def test_large_deletion1():
    name = 'large_deletion1'
    # This should just be GAG now
    query = utils.mod_hxb2[123:1626]
    utils.write_fasta({'gag': query}, cwd / f'{name}.fasta')
    # I also verified this sequence by using LANL sequence locator
    gag = valid_genes['gag']
    # Just sanity check that my query sequence is GAG
    assert gag == query
    # Align the query to the target
    # I have to skip this section because minimap2 does not run on windows
    # Otherwise I could call it as a subprocess
    # For now I will manually align the query to the target
    alignment = utils.load_samfile(cwd / f'{name}.sam')
    # Get the genes
    # results = utils.splice_genes(query, utils.mod_hxb2, alignment,
    #                              utils.mod_annot)
    results, sequences = utils.splice_aligned_genes(query, utils.mod_hxb2,
                                                    alignment, utils.mod_annot)
    # import pdb; pdb.set_trace()
    genes = utils.coords_to_genes(results, query)
    for gene in valid_genes:
        try:
            assert genes[gene] == valid_genes[gene]
        except AssertionError as e:
            print(f'Gene {gene} does not match')
            raise AssertionError(e)
