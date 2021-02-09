import gene_splicer.gene_splicer as gene_splicer
import gene_splicer.primer_finder as primer_finder
import gene_splicer.utils as utils


def main():
    data = primer_finder.main()
    for _file in data['fasta_files']:
        gene_splicer.run(_file, outdir=data['args'].outpath)
    utils.generate_table_precursor(name=data['args'].name,
                                   outpath=data['args'].outpath)


if __name__ == '__main__':
    main()
