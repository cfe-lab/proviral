import gene_splicer
import primer_finder


def main():
    fasta_files = primer_finder.main()
    for _file in fasta_files:
        gene_splicer.run(_file)


if __name__ == '__main__':
    main()