import logging

import cfeproviral.gene_splicer as gene_splicer
import cfeproviral.primer_finder as primer_finder
import cfeproviral.utils as utils


def main():
    logging.basicConfig(level=logging.WARNING)
    data = primer_finder.main()
    for _file in data['fasta_files']:
        gene_splicer.run(_file, outdir=data['args'].outpath)
    utils.generate_table_precursor(name=data['args'].name,
                                   outpath=data['args'].outpath)


if __name__ == '__main__':
    main()
