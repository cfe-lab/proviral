import logging
import re
import Levenshtein
import gene_splicer.utils as utils
from gotoh import align_it

logger = logging.getLogger(__name__)


class PrimerFinder:
    def __init__(self,
                 full_sample,
                 primer,
                 direction,
                 hxb2_start,
                 hxb2_end,
                 validation_size=25,
                 sample_size=50,
                 min_primer_length=2) -> None:
        self.primer = primer
        self.primer_length = len(self.primer)
        self.full_sample = full_sample
        self.sample_size = sample_size
        self.direction = direction
        self.get_sample()
        self.hxb2_start = hxb2_start
        self.hxb2_end = hxb2_end
        self.start = None
        self.sample_primer = None
        self.is_valid = False
        self.aln = None
        self.validation_size = validation_size
        self.min_primer_length = min_primer_length
        self.is_full_length = False
        self.find_longest_primer()
        self.validate()

    def get_sample(self):
        if self.direction == 'fwd':
            self.sample = self.full_sample[:self.sample_size]
        else:
            self.sample = self.full_sample[-self.sample_size:]

    @staticmethod
    def expand_mixtures(seq):
        old_mixtures = {''}
        for mixture in seq:
            new_mixtures = set()
            for nuc in utils.mixture_dict.get(mixture, mixture):
                for old_mixture in old_mixtures:
                    new_mixtures.add(old_mixture + nuc)
            old_mixtures = new_mixtures
        return old_mixtures

    def find_longest_primer(self):
        if self.direction == 'fwd':
            hxb2_start_offset = -1
        else:
            hxb2_start_offset = 0
        hxb2_end_offset = 0
        # Minus min_primer_length because the smallest primer we want to try
        # looking for is min_primer_length nucleotides
        for i in range(len(self.primer) - self.min_primer_length):
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

    def validate(self):
        if not self.sample_primer:
            return
        elif len(self.sample_primer) >= 7:
            self.is_valid = True
            if len(self.sample_primer) >= self.primer_length - 3:
                self.is_full_length = True
            return
        else:
            sample_slice, hxb2_slice = self.get_slices()
            if not sample_slice:
                return
            self.aln = self.align(hxb2_slice, sample_slice)
            logger.debug('\n' + self.aln['aligned_target'] + '\n' +
                         self.aln['aligned_query'])
            self.validate_alignment()
            if self.aln['is_valid']:
                self.is_valid = True
                return
        return

    def validate_alignment(self):
        if self.direction == 'fwd':
            # Take the total length and subtract the length not including leading dashes, result will be how many leading dashes
            sample_start_gaps = len(self.aln['aligned_target']) - len(
                self.aln['aligned_target'].lstrip('-'))
            hxb2_start_gaps = len(self.aln['aligned_query']) - len(
                self.aln['aligned_query'].lstrip('-'))
            if (sample_start_gaps == 0) and (hxb2_start_gaps == 0):
                self.aln['is_valid'] = True
        else:
            sample_end_gaps = len(self.aln['aligned_target']) - len(
                self.aln['aligned_target'].rstrip('-'))
            hxb2_end_gaps = len(self.aln['aligned_query']) - len(
                self.aln['aligned_query'].rstrip('-'))
            if (sample_end_gaps == 0) and (hxb2_end_gaps == 0):
                self.aln['is_valid'] = True

    def get_slices(self):
        sample_slice = ''
        hxb2_slice = ''
        # The slice should include the primer
        if self.direction == 'fwd':
            sample_slice = self.full_sample[self.start:self.end +
                                            self.validation_size]
            hxb2_slice = utils.hxb2[self.hxb2_start:self.hxb2_end +
                                    self.validation_size]
        elif self.direction == 'rev':
            # Sometimes this will fall below 0, the max is to prevent Python from accessing negative indicies from the array
            sample_slice = self.full_sample[max(
                0,
                len(self.full_sample) - self.sample_size + self.start -
                self.validation_size):len(self.full_sample) -
                                            self.sample_size + self.end]
            hxb2_slice = utils.hxb2[self.hxb2_start -
                                    self.validation_size:self.hxb2_end]
        if len(sample_slice) == 0:
            logger.debug(
                'Sample slice size is 0! \n'
                f'start: {self.start} \n'
                f'end: {self.end} \n'
                f'hxb2_start: {self.hxb2_start} \n'
                f'hxb2_end: {self.hxb2_end} \n'
                f'full_sample: {self.full_sample} \n'
            )
            return None, None
        return sample_slice, hxb2_slice

    @staticmethod
    def align(query_seq: str, target_seq: str):
        result = {
            'score': None,
            'aligned_query': None,
            'aligned_target': None,
            'start': None,
            'is_valid': False,
            'dist': None,
            'end_dist': None,
        }
        gap_open_penalty = 15
        gap_extend_penalty = 3
        use_terminal_gap_penalty = 1
        best_acontig = best_atarget = best_target = best_score = None
        for target_nucs in PrimerFinder.expand_mixtures(target_seq):
            aligned_query, aligned_target, score = align_it(
                query_seq, target_nucs, gap_open_penalty, gap_extend_penalty,
                use_terminal_gap_penalty)
            if best_score is None or score > best_score:
                best_acontig = aligned_query
                best_atarget = aligned_target
                best_target = target_nucs
                best_score = score
        if not best_acontig:
            result['valid'] = False
            return None
        aligned_query = best_acontig
        aligned_target = best_atarget
        target_nucs = best_target
        result['score'] = best_score
        match = re.match('-*([^-](.*[^-])?)', aligned_target)
        result['aligned_query'] = aligned_query
        result['aligned_target'] = aligned_target
        result['start'] = match.start(1)
        end = match.end(1)
        result['query_match'] = aligned_query[result['start']:end].replace(
            '-', '')
        result['dist'] = Levenshtein.distance(target_nucs,
                                              result['query_match'])
        stripped_contig = aligned_query.lstrip('-')
        overhang = len(aligned_query) - len(stripped_contig)
        if overhang > 0:
            stripped_target = target_nucs[overhang:]
            result['end_dist'] = Levenshtein.distance(stripped_target,
                                                      result['query_match'])
        else:
            stripped_contig = aligned_query.rstrip('-')
            overhang = len(aligned_query) - len(stripped_contig)
            if overhang == 0:
                result['end_dist'] = result['dist']
            else:
                stripped_target = target_nucs[:-overhang]
                result['end_dist'] = Levenshtein.distance(
                    stripped_target, result['query_match'])
        return result

    def __str__(self):
        return '\t'.join((str(x) for x in (self.start, self.direction,
                                           self.sample_primer, self.is_valid)))

    def __repr__(self):
        return '\t'.join((str(x) for x in (self.start, self.direction,
                                           self.sample_primer, self.is_valid)))
