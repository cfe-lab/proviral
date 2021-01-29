from gene_splicer import primer_finder
import gene_splicer.gene_splicer as gene_splicer
import gene_splicer.primer_finder as primer_finder


def main():
    data = primer_finder.main()
    for _file in data['fasta_files']:
        gene_splicer.run(_file, args=data['args'])
    if not data['fasta_files']:
        data['args'].table_precursor_csv.touch()
        data['args'].aligned_table_precursor_csv.touch()


if __name__ == '__main__':
    main()
