mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}

def reverse_and_complement(seq):
    return ''.join(complement_dict[nuc] for nuc in reversed(seq))

def read_fasta(fasta_file):
    with open(fasta_file) as fasta:
        name, seq = None, []
        for line in fasta:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            yield (name, ''.join(seq))


hxb2 = next(read_fasta(cwd / 'hxb2.fasta'))[1]