import pdb
import utils


class PrimerFinder:
    def __init__(self,
                 primer,
                 sample,
                 direction,
                 hxb2_start,
                 hxb2_end,
                 validation_size=6) -> None:
        self.primer = primer
        self.primer_length = len(self.primer)
        self.sample = sample
        self.direction = direction
        self.hxb2_start = hxb2_start
        self.hxb2_end = hxb2_end
        self.start = None
        self.sample_primer = None
        self.is_valid = False
        self.find_longest_primer()
        self.validate(validation_size)

    def expand_mixtures(self, seq):
        old_mixtures = {''}
        for mixture in seq:
            new_mixtures = set()
            for nuc in utils.mixture_dict.get(mixture, mixture):
                for old_mixture in old_mixtures:
                    new_mixtures.add(old_mixture + nuc)
            old_mixtures = new_mixtures
        return old_mixtures

    def find_longest_primer(self):
        # Minus 2 because the smallest primer we want to try looking for is 2 nucleotides
        hxb2_start_offset = 0
        hxb2_end_offset = 0
        for i in range(len(self.primer) - 2):
            # print('i', i)
            primer_substring = None
            if self.direction == 'fwd':
                primer_substring = self.primer[i:]
                hxb2_start_offset += 1
            elif self.direction == 'rev':
                if i == 0:
                    primer_substring = self.primer[:]
                else:
                    primer_substring = self.primer[:-i]
                hxb2_end_offset -= 1
            else:
                raise ValueError('direction must be either "rev" or "fwd"')

            # print('primer_substring', primer_substring)
            primer_substrings = self.expand_mixtures(primer_substring)
            for sub in primer_substrings:
                try:
                    self.start = self.sample.index(sub)
                    self.hxb2_start += hxb2_start_offset
                    self.hxb2_end += hxb2_end_offset
                    self.end = self.start + len(sub)
                    self.sample_primer = sub
                    return
                except ValueError:
                    continue

    def validate(self, validation_size):
        if not self.sample_primer:
            return
        sample_slice = ''
        hxb2_slice = ''
        if self.direction == 'fwd':
            sample_slice = self.sample[self.end:self.end + validation_size]
            hxb2_slice = utils.hxb2[self.hxb2_end:self.hxb2_end +
                                    validation_size]
        elif self.direction == 'rev':
            sample_slice = self.sample[self.start - validation_size:self.start]
            hxb2_slice = utils.hxb2[self.hxb2_start -
                                    validation_size:self.hxb2_start]
        sample_slice_len = len(sample_slice)
        if sample_slice_len == 0:
            return
        self.distance = 0
        for i, nuc in enumerate(sample_slice):
            if nuc != hxb2_slice[i]:
                self.distance += 1
        percentage_match = (sample_slice_len -
                            self.distance) / sample_slice_len
        if (self.distance <= 1) or (percentage_match > 0.75) or (len(
                self.sample_primer) > 6):
            self.is_valid = True

    def __str__(self):
        return '\t'.join(
            (str(x) for x in (self.start, self.direction, self.sample_primer)))

    def __repr__(self):
        return '\t'.join(
            (str(x) for x in (self.start, self.direction, self.sample_primer)))