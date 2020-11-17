import enum
import pdb
import pytest
import math
import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

import gene_splicer.utils as utils
import Levenshtein
from gene_splicer.probe_finder import ProbeFinder
from gene_splicer.probe_finder import TARGET_SEQUENCES
from gene_splicer.primer_finder import primers, mixture_dict
from gene_splicer.utils import hxb2


@pytest.fixture
def data():
    return {'cwd': Path(os.path.realpath(__file__)).parent}


# def test_3prime_end(data):
#     # finder = ProbeFinder(data['seq'], TARGET_SEQUENCES['round2_rev_primer'])
#     hxb2_target_start = primers['rev']['hxb2_start'] - 100
#     hxb2_target_end = primers['rev']['hxb2_end']
#     hxb2_target_seq = utils.hxb2[hxb2_target_start:hxb2_target_end]
#     finder = ProbeFinder(hxb2_target_seq, data['seq'])
#     finder.start += hxb2_target_start
#     print()
#     print(finder)
#     validated = validate_primer(finder, data['seq'], primers['rev'])
#     print()
#     print(validated)


def validate_primer(finder, finder_seq, target, tolerance=1):
    if finder.is_reversed:
        finder_seq = utils.reverse_and_complement(finder_seq)
    error = None
    primer_in_finder = None
    finder_seqsize = len(finder_seq)
    matched_finder_size = len(finder.contig_match)
    overhang = matched_finder_size - finder_seqsize
    primer_in_finder_hxb2_start_coord = max(finder.start, target['hxb2_start'])
    primer_in_finder_hxb2_end_coord = min(finder.start + finder_seqsize,
                                          target['hxb2_end'])
    # Get the coordinates of the primer relative to the finder sequence
    primer_in_finder_start_coord = primer_in_finder_hxb2_start_coord - finder.start
    primer_in_finder_end_coord = primer_in_finder_hxb2_end_coord - finder.start
    if ((target['hxb2_start'] == 9603) and
        (primer_in_finder_end_coord >= matched_finder_size - overhang)):
        primer_in_finder_start_coord -= overhang
        primer_in_finder_end_coord -= overhang

    # Get the primer sequence of the finder sequence
    primer_in_finder = finder_seq[
        primer_in_finder_start_coord:primer_in_finder_end_coord]

    # Get the sequence of the true primer that overlaps the finder sequence
    primer_start_coord = max(
        0, primer_in_finder_hxb2_start_coord - target['hxb2_start'])
    primer_end_coord = primer_start_coord + len(primer_in_finder)
    real_primer = target['seq'][primer_start_coord:primer_end_coord]
    result = {
        'finder_seq': finder_seq,
        'target_seq': primer_in_finder,
        'contig_hxb2_start': finder.start,
        'contig_hxb2_end': finder.start + matched_finder_size,
        'seq_start': primer_in_finder_start_coord,
        'seq_end': primer_in_finder_end_coord,
        'hxb2_start': primer_in_finder_hxb2_start_coord,
        'hxb2_end': primer_in_finder_hxb2_end_coord,
        'overhang': overhang,
        'contig_match': finder.contig_match,
        'dist': 0,
        'finder_dist': finder.dist,
        'real_primer': real_primer,
        'error': error
    }

    if len(real_primer) != len(primer_in_finder):
        result[
            'error'] = 'Real primer did not match length of primer in finder'
        return result
    if not real_primer:
        result['error'] = 'real primer not found at expected coordinates'
        return result
    elif not primer_in_finder:
        result[
            'error'] = 'primer in contig sequence not found at expected coordinates'
        return result
    for i in range(len(real_primer)):
        try:
            our_nuc = primer_in_finder[i]
        except IndexError:
            pass
        real_nuc = real_primer[i]
        if our_nuc != real_nuc:
            if real_nuc in utils.mixture_dict and our_nuc in utils.mixture_dict[
                    real_nuc]:
                pass
            else:
                result['dist'] += 1
        if result['dist'] > tolerance:
            result['error'] = 'mismatches in primer > tolerance'
            return result
    # if finder.dist > (finder_seqsize / 10):
    #     result['error'] = f'Poor alignment (finder_dist > {finder_seqsize/10})'
    #     return result
    return result


def primer_starts_after_seq_end(finder, hxb2_segment_start):
    if hxb2_segment_start > finder.start:
        return True
    else:
        return False


def primer_ends_before_seq_start(finder, hxb2_segment_end):
    if finder.start > hxb2_segment_end:
        return True
    else:
        return False


def primer_overlaps_seq(finder, hxb2_segment_start, hxb2_segment_end):
    if hxb2_segment_start <= finder.start <= hxb2_segment_end:
        return True
    else:
        return False


def get_primer_from_sample(sample_probe, primer_buffer, hxb2_segment,
                           hxb2_segment_start, hxb2_segment_end, primer):
    """Consider an segment of hxb2 H of length L(H) containing a primer P. The finder is H aligned to sample_probe. Attempt to locate a primer within sample_probe that overlaps any portion of P.

    Args:
        sample_probe (string): A pre-defined segment of the sample sequence
        hxb2_segment (string): A piece of HXB2 containing either the second round fwd primer or rev primer
        hxb2_segment_start (int): The start position (0-based) in hxb2 of the segment
        primer (dict): A dictionary containing primer details. e.g.
            {
                'seq': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
                'nomix': 'TAAGCCTCAATAAAGCTTGCCTTGAGTGC',
                'hxb2_start': 9604 - 1,
                'hxb2_end': 9632,
                'direction': 'rev'
            }
        tolerance (int, optional): The number of mismatches to accept between the primer located within the sample and the true primer. Defaults to 1.
    """
    # Align hxb2_segment containing primer to seq
    alignment = ProbeFinder(sample_probe, hxb2_segment)
    print(alignment.aligned_target + '\n' + alignment.aligned_contig)
    # Check if the alignment was valid and if not return an error
    if not alignment.valid:
        return {'error': 'Too many sequences to unpack'}

    sample_probe_primer_start = alignment.start
    if primer['direction'] == 'rev':
        sample_probe_primer_start += primer_buffer
    sample_probe_primer_end = sample_probe_primer_start + len(primer['seq'])
    sample_probe_primer_seq = alignment.aligned_contig[
        sample_probe_primer_start:sample_probe_primer_end]
    return sample_probe_primer_seq
    # Adjust for hxb2 coordinates
    # alignment_hxb2_start = alignment.start + hxb2_segment_start
    # if primer_overlaps_seq(alignment, hxb2_segment_start, hxb2_segment_end):
    #     result['overlaps'] = True
    # elif primer_starts_after_seq_end(alignment, hxb2_segment_start):
    #     result['error'] = f'{end} sample probe ends before hxb2 primer start'
    # elif primer_ends_before_seq_start(alignment, hxb2_segment_end):
    #     result['error'] = f'{end} sample probe starts after hxb2 primer end'

    # Get the start coordinate of the primer within the sample segment relative to hxb2
    # primer_in_sample_probe_hxb2_start_coord = max(alignment_hxb2_start,
    #                                               primer['hxb2_start'])
    # primer_in_sample_probe_hxb2_end_coord = min(
    #     alignment_hxb2_start + probelen, primer['hxb2_end'])

    #
    # print()
    # print('sample_probe', sample_probe)
    # print('hxb2_segment', hxb2_segment)
    # print('primer', primer)
    # print('alignment', alignment)


# def test_fwd():
#     probe_length = 50
#     primer_buffer = 29
#     fwd_sample_probe = 'GCGCCCGAACAGGGACCTGAAAGCGAAAGCCCCCCCCCCCCCCCCCCCCC'
#     fwd_hxb2_start = primers['fwd']['hxb2_start']
#     fwd_hxb2_end = primers['fwd']['hxb2_end'] + primer_buffer
#     fwd_hxb2_seq = hxb2[fwd_hxb2_start:fwd_hxb2_end]
#     fwd_primer = get_primer_from_sample(fwd_sample_probe, primer_buffer,
#                                         fwd_hxb2_seq, fwd_hxb2_start,
#                                         fwd_hxb2_end, primers['fwd'])

# def test_mismatch_revprimer(data):
#     sample_probe = None
#     direction = 'rev'
#     primer_buffer = 29
#     probelen = 100
#     with open(data['cwd'] / 'data' / 'real_example.txt') as f:
#         sample_probe = f.read().strip()[-100:]
#     hxb2_target_start = primers[direction]['hxb2_start'] - primer_buffer
#     hxb2_target_end = primers[direction]['hxb2_end']
#     hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]

# def test_mismatch_revprimer(data):
#     seq = None
#     name = 'rev'
#     primerbuffer = 20
#     # probelen = 100
#     # with open(data['cwd'] / 'data' / 'mismatched_revprimer.txt') as f:
#     with open(data['cwd'] / 'data' / 'real_example.txt') as f:
#         seq = f.read().strip()
#     hxb2_target_start = primers[name]['hxb2_start'] - primerbuffer
#     hxb2_target_end = primers[name]['hxb2_end']
#     hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
#     print(hxb2_target_seq)
#     print(data)
#     # truncate seq to appropriate size
#     # seq = seq[-probelen:]

#     # finder = ProbeFinder(hxb2_target_seq, seq)
#     finder = ProbeFinder(seq, hxb2_target_seq)
#     print(finder)
#     finder.start += hxb2_target_start
#     print(finder)
#     validated = validate_primer(finder, seq, primers[name])
#     print(validated)


def generate_test_cases():
    for fwd_size in range(0, 29):
        for rev_size in range(0, 29):
            for fwd_padding in range(0, 30, 10):
                for rev_padding in range(0, 30, 10):
                    for fill_len in range(100, 1000, 200):
                        base = 'A'
                        fwd_pad = base * fwd_padding
                        rev_pad = base * rev_padding
                        fwd_primer = primers['fwd']['nomix'][-fwd_size:]
                        if fwd_size == 0:
                            fwd_primer = ''
                        rev_primer = primers['rev']['nomix'][:rev_size]
                        if rev_size == 0:
                            rev_primer = ''
                        filler = utils.hxb2[
                            primers['fwd']['hxb2_end']:
                            primers['fwd']['hxb2_end'] +
                            math.ceil(fill_len / 2)] + utils.hxb2[
                                primers['rev']['hxb2_start'] -
                                math.ceil(fill_len /
                                          2):primers['rev']['hxb2_start']]
                        sequence = fwd_pad + fwd_primer + filler + rev_primer + rev_pad
                        expected = True
                        if fwd_size < 4 or rev_size < 4 or fwd_padding > 29 or rev_padding > 29:
                            expected = False
                        yield {
                            'fill_len': fill_len,
                            'fwd_padding': fwd_padding,
                            'fwd_size': fwd_size,
                            'rev_size': rev_size,
                            'rev_padding': rev_padding,
                            'sequence': sequence,
                            'expected': expected
                        }


def compute_primer_distance(seq, primer_seq, min_primer=2):
    dist = 0
    gap_count = 0
    for i, nuc in enumerate(seq):
        target_nuc = primer_seq[i]
        if nuc != target_nuc:
            if target_nuc in mixture_dict and nuc in mixture_dict[target_nuc]:
                pass
            elif nuc == '-':
                gap_count += 1
            else:
                dist += 1
    if len(seq) - gap_count <= min_primer:
        return gap_count
    return dist


# def test_case():
#     case = {
#         'fwd_padding': 0,
#         'fwd_size': 6,
#         'rev_size': 6,
#         'rev_padding': 0,
#         'sequence': 'CGAAAGAAAAAAAAAAAAAAAAAAAATAAGCC',
#         'expected': True
#     }
#     probe_length = 50
#     primer_buffer = 29
#     fwd_sample_probe = case['sequence'][:probe_length]
#     rev_sample_probe = case['sequence'][-probe_length:]
#     fwd_hxb2_start = primers['fwd']['hxb2_start']
#     fwd_hxb2_end = primers['fwd']['hxb2_end'] + primer_buffer
#     fwd_hxb2_seq = hxb2[fwd_hxb2_start:fwd_hxb2_end]
#     rev_hxb2_start = primers['rev']['hxb2_start'] - primer_buffer
#     rev_hxb2_end = primers['rev']['hxb2_end']
#     rev_hxb2_seq = hxb2[rev_hxb2_start:rev_hxb2_end]
#     import pdb
#     pdb.set_trace()
#     fwd_primer = get_primer_from_sample(fwd_sample_probe, primer_buffer,
#                                         fwd_hxb2_seq, fwd_hxb2_start,
#                                         fwd_hxb2_end, primers['fwd'])
#     rev_primer = get_primer_from_sample(rev_sample_probe, primer_buffer,
#                                         rev_hxb2_seq, rev_hxb2_start,
#                                         rev_hxb2_end, primers['rev'])
#     print()
#     print(case)
#     print(f'fwd_primer: {fwd_primer.sample_primer}')
#     print(f'rev_primer: {rev_primer.sample_primer}')
#     if fwd_primer.sample_primer and rev_primer.sample_primer:
#         print('matched expectations?', True == case['expected'])
#     else:
#         print('matched expectations?', False == case['expected'])


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
            for nuc in mixture_dict.get(mixture, mixture):
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
                    # self.hxb2_start -= 1
                    self.hxb2_start += hxb2_start_offset
                    self.hxb2_end += hxb2_end_offset
                    self.end = self.start + len(sub)
                    self.sample_primer = sub
                    return
                except ValueError:
                    continue

    def validate(self, validation_size):
        sample_slice = ''
        hxb2_slice = ''
        if self.direction == 'fwd':
            sample_slice = self.sample[self.end + 1:self.end + 1 +
                                       validation_size]
            hxb2_slice = utils.hxb2[self.hxb2_end + 1:self.hxb2_end + 1 +
                                    validation_size]
        elif self.direction == 'rev':
            sample_slice = self.sample[self.start - 1 -
                                       validation_size:self.start - 1]
            hxb2_slice = utils.hxb2[self.hxb2_start - 1 -
                                    validation_size:self.hxb2_start - 1]
        sample_slice_len = len(sample_slice)
        if sample_slice_len == 0:
            return
        self.distance = 0
        for i, nuc in enumerate(sample_slice):
            if nuc != hxb2_slice[i]:
                self.distance += 1
        percentage_match = (sample_slice_len -
                            self.distance) / sample_slice_len
        import pdb
        pdb.set_trace()
        if (self.distance <= 1) or (percentage_match > 0.75):
            self.is_valid = True

    def __str__(self):
        return '\t'.join(
            (str(x) for x in (self.start, self.direction, self.sample_primer)))

    def __repr__(self):
        return '\t'.join(
            (str(x) for x in (self.start, self.direction, self.sample_primer)))


def test_case():
    case = {
        'fill_len': 100,
        'fwd_padding': 0,
        'fwd_size': 4,
        'rev_size': 4,
        'rev_padding': 0,
        'sequence':
        'AAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAG',
        'expected': True
    }
    case = {
        'fill_len': 100,
        'fwd_padding': 10,
        'fwd_size': 1,
        'rev_size': 4,
        'rev_padding': 0,
        'sequence':
        'AAAAAAAAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAG',
        'expected': False
    }
    case = {
        'sequence':
        'GGGACCTGAAAGCGAAAGAAAAACCAAAGAAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACAGCAAGAGGCGAGGGGCGGCGACTGGTAAGTACGCCAAAAATTTTAACTAGCAAAAGCTAAAAGAAAAGAAATGAGTGCAAAAGCGTCAATATTAAGTGGGGGAAAATTAAATAAATGGGAAAAAATTCGGTTAAGGCCAGGGGAAAAAAAAAAGTATAAATTAAAACACATCGTATGGGCAAGCAGGAAGCTAAAACAATTCGCAGTTAACCCTGGCCTGTTAAAAACAGCAAAAGGCTGTAGACAAATACTGGGACAATTACAACCGTCCCTTCAGACAGAATCAGAAAAACTTAAATCCTTATACAATACAGTAGCAACACTCTATTGTGTGCATCAAGAGATAGAAATAAAAAACACCAAGAAAGCTTTAGACAAAATAAAGAAAAAACAAAACAAAAGTAAGAAAAGGGCACAGCAAGCAGCAGCTAACACAGGAAACAAAGGCCAGGTCAGCCAAAATTTCCCTATAGTACAAAACCTCCAGGGACAAATGGTACATCAGGCTATATCACCTAAAACTTTAAATGCATGGGTAAAAGTAATAAAAAAGAAAGCTTTTAGCCCAGAAGTAATACCCATGTTTACAGCACTATCAGAAAAAGCCACTCCACAAAATTTAAACACCATGCTAAACACAGTAGAGAAACATCAAGCAGCCATGCAAATGTTAAAAAAAACCATCAATAAAAAAGCTGCAGAATGAAATAGATTGCATCCAGTGCAGGCAGAACCCGTTGCACCAGGTCAAATAAAAAACCCAAGGGGAAGTGACATAGCAGAAACTACTAGTACCCTTCAGAAACAAATAGCATAAATAACACATAATCCACCTATCCCAGTAGAAGACATCTATAAAAAATAAATAATCCTGAGTTTAAATAAAATAGTAAAAATGTATAGTCCTGCCAGCATTCTAGACATAAGACAAGGGCCAAAGAAACCTTTTAAAGACTATGTAAACCGGTTCTATAAAACTCTAAAAGCCAAGCAAGCCACACAAAAGGTAAAAAATTAAATAACAGAAACCTTGTTAGTCCAGAATGCAAACCCAAATTGTAAAACTATTTTAAAAGCATTAGGTCCAGAAGCTACACTAAAAAAAATAATAACAGCATGTCAAGAAGTAAAAAAACCCAACCACAAAGCAAGGGTCTTGGCAGAAGCAATAAGCCAAGCAACAACTAATGCAGCCGTAATAATGCAAAAAGGCAATTTTAAAAACCAAAAAAAAAATTGTTAAATGTTTCAATTGTGGCAAAAAAAGGCATATAGCTAAAAATTGCAAGGCCCCTAAAAAAAAGAGGTGCTAAAAATGTAAAAAAAAAAAACACCAAATAAAAAACTGTACTAAAAAACAGGCTAATTTTTTAAGAAAAATCTAGCCTTCCCGCAAGGGGGGGCCAGGAAATTTCATTCAGAGCAAACTAAAGCCAACAGCCCCACCAAAAAAAACCTTCAGGTTTAGAAAAGAAACAACAACTCCCTCTCAGAAGCAGAAAACAACAAACAAGAACAACATGTATCCCTTAGCCTCCCTCAAATCACTCTTTAGCAACAACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAAAAAGCTCTATTAAATACAGAAGCAAATAATACCGTATTAAAAAATATAACTTTGCCAGAAAAATAAAAACCAAAAATAATAAAAAAAATTAAAGGTTTTATCAAAGTAAAACAGTATAATCAAATACCTATAAAAATCTGTAAACATAAAGCTATAAGTACAGTGTTAGTAAAACCCACACCTGTCAACATAATTAAAAAAAATCTGTTAACTCAAATTAGCTGCACTTTAAATTTCCCCATTAGTCCTATAAACACTGTACCAGTAAAATTAAAGCCAGAAATAAATAGCCCAAAAGTTAAACAATAGCCATTAACAGAAAAAAAAATAAAAGCATTAATAAAAATCTGTACAAAAATAAAAAAAAAAAAAAAATTTCAAAAATTAAGCCTAAAAATCCATACAATACTCCAATATTTGCCATAAAAAAAAAAAACGGTACCAAATAAAAAAAATTAGTAAATTTTAAAAAACTTAATAAAAAAACTCAAAACTTCTAAAAAGTTCAATTAAAAATACCACATCCCGCAAGGTTAAAAAAAAGTAAGTCAGTAACGGTGCTAAACGTAAGTAATGCATATTTTTCAGTTCCCTTAAATAAAAATTTCAGAAAGTATACTGCATTCACCATACCTAATATAAACAATAAAACACCAAGAATTAAATACCAGTACAATGTGCTTCCACAGGAATAAAAAAAATCACCAGCAATATTCCAAAGTAGCATAACAAAAATCTTAAAGCCTTTTAAAAAACAAAACCCAAACATAACTATCTATCAATACATAAATAATCTGTATGTAGCATCTAACTTAAAAATAAAGCAACATAAAACAAAAATACAGAAACTAAAAAATCATCTGTTAAAATAAAAGTTTCTTACACCAGACAAAAAACATCAAAAAAAACCTCCATTCCTTTAAATAGAGTATAAACTCCATCCTAATAAATAAACAGTACAGCCTATAGTGCTGCCAGAAAAAAACAGCTAAACTGTCAATAATATACAGAAGTTAGTAAAAAAATTAAATTAGGCAAGTCAAATTTACCCAGAAATTAAAGTAAGGCAATTATGTAAACTCCTTAAAAAAACCAAAGCACTAACAAAAGTAGTACCACTAACAGCAAAAGCAAAGCTAAAACTGGCAAAAAACAGGAAAATTCTAAAAAAACCAGTACATAAAGTGTACTATAACTCATCAAAAAAATTAATAGCAAAAATACAAAAGCAGGAGTTAAACCAATAGACATATCAAATCTATCAAAAGCCAGGTAAAAATCTAAAAACAGAAAAATATGCAAAAATAAAAAGTGCCCACACTAATAATGTAAAACAATTAACAAAGGCAGTGCAAAAAATAGCCACAAAAAGCATAGTAATATAAAAAAAAACTCCTAAATTTAAACTACCCATACAAAAAAAAACATAGAAAACATAGTAAACAAAGTATTAGCAGGCCACCTAAATTCCTAAATAAAAATTTGTCAATACCCCTCCCTTAGTAAAATTATAGTACCAGTTAAAAAAAAAACCCATAAAAAAAGCAGAAACCTTCTATGTAAATAAGGCAGCTAACAAGAAAACTAAAGCAGAAAAAGCAGAATATGTTACTAACAAAAAAAAACAAAAGATTATCTCCATAACTAACACAACAAATCAAAAAACTAAGTTACAAGCAATTCATCTAGCTTTGCAGAATTCAGAACCAAAAATAAACATAGTAACAAACTCACAATATGCATTAAAAATCATTCAAGCACAACCAAATAAAAGTAAATCAAAGTTAGTCAGTCAAATAATAAAAAAGTTAATAAAAAAAAAAAAAATCTATTTAGCATAAATGCCAGCACATAAAAAAATTAAAAAAAATAAACAAGTTAATAAATTAGTCAGTGCTAAAATCAGAAAAGTACTATTTTTAAATAAAATAAACAAAGCCCAAAAAAACCATAAAAAATATCATAGCAATTAAAAAGCTATAGTTAATAATTTTAACCTGCCACCTGTAGTAGCAAAAAAAATAGTAGCCAGCTGTAATAAATGTCAGCTAAAAAAAAAAGCCATACATGGCCAAGTAAACTGTAGTCCAGAAATGTAGCAACTAAATTGTACGCATTTAAAAAAAAAAGTTATCCTAGTAGCAGTTCATGTAGCCAGTAAATACATAAAAGCAGAAGTTATCCCAGCAAAAACAGAACAAAAAACAGCATACTTTCTCTTAAAATTAGCAAAAAAATAGCCAGTAAAAACAATACATACAAACAATAGCCCCAATTTTACCAGCACTGCAGTTAAAGCCGCCTATTAGTAAGCAGAAATCAAGCAAAAATTTAAAGTTCCCTACAATCCCCAAAGTCAAAAAGTAGTAAAATCTATAAATAAAAAATTAAAAAAAATTATAAAACAGGTAAAAAATCAAGCTAAGCATCTTAAAACAGCAGTACAAATAGCAGTACTCATCCACAATTTTAAGAGAAAAGGGGGGATTGGGGGATACACAGCAGGGAAAAGAATAATAGACATAATAGCAACAAACATACAAACTACAGAATTACAAAAACAAATTACAAAACTTCAAAATTTTCAGGTTTATTACAGAGACAACAAAGATCCACTTTGGAAAGAACCAGCCAAGCTTCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAACATAATAGTGACATAAAAGTAGTGCCAAGAAAAAAAGCAAAGATCATTAGGCATTATGAAAAACAAATGGCAGGCAATAATTGTGTGGCAAGTAGACAGGATAAGGATTAAAACATGAATAAGTTTAGTAAAACACCATATGTATATTTCAAAGAAAGCTAAAAGATGGGTTTATAAACATCACTATGAAAACCCTCATCCAAGAGTAGGTTCAAAAGTACACATCCCACTAGGGGATGCTAAATTGGTAGTAACAACATATTGAGGTCTGCACACAGGAGAAAAAAACTGGCATTTGGGTCAGGGAGTCTCCATAAAATGAAGGAAGAGGGAATATAGCACGCAAGTAGACCCTGGCCTAGCAGACCAACTAATTCACCTGTATTACTTTAATTGTTTTTCAAACTCTGCTATAAAAAAGGCCATATTAGGACATATAGTTAGTCCTAGTTGTAAATATCCAACAGAACATAACAAGGTAGAATCTCTACAGTACTTGGCACTAACAGCATTAATAACACCAAAAAAAATAAAGCCACCTTTGCCTAGTGTTATAAAACTGACAGAGGATAGATGGAACGAGCCCCAGAAGACCAAGGGCCACCAAAGAAGCCATACCATAAATGGGCACTAGAGCTTTTAAAGAAACTTAAAAATAAAGCTGTCAGACATTTTCCCAGAATGTGGCTCCATGGCTTAGGGCAATACATCTATGAAACCTATAAGAATACTTGGGCAGGAGTGAAAGCCATAATAAGAATTCTGCAACAACTACTGTTTATTCATTTCAGAATTAGGTGTCGACATAGCAGAATAGGCATTACTCTACAGAGGAGAACAAGAAATAAAGCCAGTAAATCCTAAACTAGAGCCTTAGAAGCATCCAGGAAGTCAGCCTAAGACTGCTTGTACCCCTTGCTATTGTAAAAAGTGTTGCTTTCACTGTCAAGTTTGTTTCATAAAAAAAGGCTTAGGCATCTCCTATGGCAGGAAAAAGCGGAGACAGCGACAAAGAACTCCTCAGGACAATCAGGCTCATCAAGTTCCTCTACCAAAGCAGTAAGTAATATATGTAATGCAATCTTTAAAAATAGTATCAATAGCAGCCTTAGTAGTAGCAGCAATAATAGCAATAGTTGTGTAAACCATAGTATACATAAAATATAAAAAAATATTAAGACAAAAAAAAATAAATAAGCTAATTAAAAAAATAAGTGAAAAAGCAAAAAACAGTGGCAATAAAAGTAAAGAAAATCAGAAGAAATTATCAGCACTTGTGAAAATGGGGCACCAGGCTCACAATGCTCCTTGAAATGTTAATAATCTGTAATGCTACAGAACAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTAAAAAAAAGCAATCACCACTCTATTTTGTGCATCAAATGCTAAGGCATATAATACAAAGGTACATAATGTTTAAGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAAAAGTAAAATTAAAAAATGTAACAGAAAATTTTAATATGTAAAAAAATAACATGGTAAAACAGATGCATAAAAATATAATTAGTTTATAGAATCAAAGCCTAAAGCCATGTGTCAAATTAACCCCACTCTGTGTTACTCTACATTGCACTAACGCGAATAATACTGACGCGAATAATACTAGTAATAAAAATATTACTTTAAAAATAGAAAAAGAAAAAATACAAAACTGCTCTTTCAATGTCACTTCAAGCATAAGAAATAAAGTGCAAAAAAAATATGCGCTCTTTTATAAACTTAATGTAGTACCAATAAATAATAATACCAACTCTAATTACAGCTCTTATAATACCAACTCTAATTACAGCTCTTATAGTTACAGCTCTTATAGGTTAATAAATTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAAGTAACCTTTAAGCCAATTCCTATACATTATTGTGCCCCGGCTAGTTTTGCAATTCTAAAGTGTAAAAATAAAAAGTTCAATGGAACAGGGCCATGTACAAATGTTAGCACAGTACAATGTACTCATGGAATTAAGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGCCTAGCAGAAAAAAAGGTAATAATTAAATCTAAAAATTTCTCAAACAATGCTAAAACCATAATAGTGCAGTTAAAAAAATCTGTAAAAATTACTTATATAAAACCTAATAACAATACAAGAAAAAGTATACCTATAAAACCAGGAAAAGCATTTTATGCAACAGAAAACATAATAAAAAACATAAAAAAAGCACATTGTAATATTAATAAAGCAGAATAAAACAAAGCTTTACAACAAGTAGCTAAAAAATTAAGAAAACAATTTAATAGTACAACAATAATCTTTACTAACTCCTCAGAAAGAAACCCAGAAATCACAACGCACAGCTTTAATTGTAAAAAAAAATTTTTCTACTGCAATACAACAAACTTGTTTAATAGTATTTAAAATACAACACATGAGTTAAGTAGTACTTAAAATAGTACTAACAGTAACAATATCACACTCCAATGCAGAATAAAACAAATTATAAATATGTGGCAGAAAGTAAAAAAAGCCATGTATGCCCCTCCCATCAGGGGGCTAATTCAATGTACATCAAATATTACAGGGCTGCTATTAACAAAAAATGGTAGTGGTAATGGTACTAAAACTAATAAAACCTTCAAACCAGAAAAAAAAAATATAAAGGACAATTAAAAAAATAAACTATATAAGTATAAAATAGTAAAAATTAAACCATTAAGAGTAGCACCCACCAAGGCAAAGAGAAAAGTGGTGCAGAGAGAAAAAAAAGCAGTAGAAATAAAAGCTTTATTCCTTAAGTTCTTGGCAGCAGCAAAAAGCACTATGGGCGCAGCGTCAATAACGCTAACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAACAGCAAAACAATTTGCTAAGGGCTATTAAAGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATAAAGCAGCTCCAGGCAAAAGTCCTAGCTGTGAAAAAATACCTAAAAAATCAACAGCTCCTGAAAATTTAAGGTTGCTCTAAAAAACTCATCTGCACCACTGCCGTGCCTTAAAATAGTAGTTAAAGTAATAAGTCTCTAAAACAAATTTAAAATAACATAACCTAAATGCAGTAAAAAAAGAAAATTAACAATTACACAGGCTTAATATACAACCTAATTAAAAAATCGCAAAACCAACAAAAAAAAAATAAAAAGAAATTATTAAAATTAAATAAATAGGCAAGTTTGTAAACTTAGTTTAACATAACAAACTAGCTGTGGTATATAAAAATATTTATAATAATAGTAAAAGGCTTAATAAGTTTAAAAATAGTTTTTACTGTACTTTCAATAGTAAATAAAGTTAAGCAGAAATACTCACCATTATCGTTTCAAACCCGCCTCCCAGCCCAAGGGGAACCCAACAGGCCCAAAGAAATCAAAAAAAAAGGTAAAAAAAAAGACAGAGACAAATCCAAAAACTTAGCAAATAAATTATTAGCAATCATTTAGGTCAACCTGCGGAGCCTGTTCCTTTTCATCTACCACCACTTAAAAGACTTACTCTTAATTGCAACAAAAATTATAAAACTTCTGAAACGCAGGGAGTAAAAAATCCTCAAGTACTGGTAAAACCTCCTCCTGTATTAAAGTCAGAAACTAAAAAGTAGTGCTGTTAGCTTGTTCAACGCCACAGCCATAGCAGTAGCTAAAGGGACAAATAAGGTTATAAAAGTAGCACAAAAAGCTTTTAAAGCAGTTCTCCACATACCTAAAAAAATAAAACAAGGCTTAAAAAGGCTTTTACTATAAAATAAAAGGCAAGTGGTCAAAAAATAAAAAGCTTAAATAGTCTGCTGTAAAGAAAAAAATAAAACAAGCTAATCCAAAAGCAAAGCCAGCAGCAAATAAAGTAAAAGCAGTATCTCAAAACCTGAAAAAATATAAAGCAGTCTCAACTAGCAATACAGCACATACCAATGCTAATTGTGCCTAACTAAAAGCACAAAAGAATAAAAATGTAGGCTTTCCAGTCAAACCTCAAGTACCTTTAAAACCAATAACTTACAAAGCAGCAGTAAATCTTAGCCACTTTTTAAGAGAAAAGGGGGGACTGGAAGGGTTAATTTACTCCCAAAGAAGACAAGACATCCTTAATTTATGGGTCTATCACACACAAGGCTTCTTCCCTAATTGGCAGAATTACACACCAGGGCCAGGGATCAAATATCCACTAACCTTTAAATGGTGCTACAAGCTAGTACCAGTAGAACCAGAGGAGGTAGAAAAGGCTAATGAAGAAAAGAACAATGTCTTGTTACACCCTATAAGCCAGCATGAAATGAATAACCCTAAGAAAGAAGTGCTAGTGTAGAAGTTTAACAGCCGCCTGGCATTTCAACACGTAGCCAAAAAGAAACATCCTGAGTTCTACAAGAACTGCTGACATTACAAGACCTGCTGACAACTGCTAACATCAAGCTTTCTACAAGGAACTTTCCGCTGGGGACTTTCCAGGGAGGTGTGGCCTGGGCAGAACTGGGAAGTGGCAAGCCCTCAAATGCTGCATATAAGCAGCTGCTTTCTGCTTGTTCTGAGTCTCTCTTGTTAAACCAAATCCAAGCCCGAAAGCTCTCTGGCTAGCTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTT'
    }
    probe_length = 50
    primer_buffer = 29
    fwd_sample_probe = case['sequence'][:probe_length]
    rev_sample_probe = case['sequence'][-probe_length:]

    fwd_primer = PrimerFinder(primers['fwd']['seq'], fwd_sample_probe, 'fwd',
                              primers['fwd']['hxb2_start'],
                              primers['fwd']['hxb2_end'])
    rev_primer = PrimerFinder(primers['rev']['seq'], rev_sample_probe, 'rev',
                              primers['rev']['hxb2_start'],
                              primers['rev']['hxb2_end'])
    print()
    print(case)
    print(f'fwd_primer: {fwd_primer.sample_primer}')
    print(f'rev_primer: {rev_primer.sample_primer}')
    if (fwd_primer.is_valid and rev_primer.is_valid):
        print('matched expectations?', True == case['expected'])
    else:
        print('matched expectations?', False == case['expected'])


def test_test_cases():
    print()
    for case in generate_test_cases():
        probe_length = 50
        primer_buffer = 29
        fwd_sample_probe = case['sequence'][:probe_length]
        rev_sample_probe = case['sequence'][-probe_length:]

        fwd_primer = PrimerFinder(primers['fwd']['seq'], fwd_sample_probe,
                                  'fwd', primers['fwd']['hxb2_start'],
                                  primers['fwd']['hxb2_end'])
        rev_primer = PrimerFinder(primers['rev']['seq'], rev_sample_probe,
                                  'rev', primers['rev']['hxb2_start'],
                                  primers['rev']['hxb2_end'])
        print()
        print(case)
        print(f'fwd_primer: {fwd_primer.sample_primer}')
        print(f'rev_primer: {rev_primer.sample_primer}')
        if (fwd_primer.is_valid and rev_primer.is_valid):
            print('matched expectations?', True == case['expected'])
        else:
            print('matched expectations?', False == case['expected'])