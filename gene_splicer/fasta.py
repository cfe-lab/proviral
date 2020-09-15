from gene_splicer.file import File


class Fasta(File):
    def read(self):
        with self.open() as fasta:
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

    def __iter__(self):
        return self.read()