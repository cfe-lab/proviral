import gene_splicer
import primer_finder


def main():
    data = primer_finder.main()
    for _file in data['fasta_files']:
        gene_splicer.run(_file, outdir=data['args'].outpath)


if __name__ == '__main__':
    main()
