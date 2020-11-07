import pytest
import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent
sys.path.append(str(cwd.parent.parent))

import gene_splicer.utils as utils
from gene_splicer.probe_finder import ProbeFinder
from gene_splicer.probe_finder import TARGET_SEQUENCES
from gene_splicer.primer_finder import primers, validate_primer
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


def validate_primer2(sample_probe,
                     end,
                     hxb2_segment,
                     hxb2_segment_start,
                     hxb2_segment_end,
                     primer,
                     tolerance=1):
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
    probelen = len(sample_probe)
    result = {'primer_in_sample': None}
    # Align hxb2_segment containing primer to seq
    alignment = ProbeFinder(hxb2_segment, sample_probe)
    # Check if the alignment was valid and if not return an error
    if not alignment.valid:
        return {'error': 'Too many sequences to unpack'}
    # Adjust for hxb2 coordinates
    alignment.start += hxb2_segment_start
    if primer_overlaps_seq(alignment, hxb2_segment_start, hxb2_segment_end):
        result['overlaps'] = True
    elif primer_starts_after_seq_end(alignment, hxb2_segment_start):
        result['error'] = f'{end} sample probe ends before hxb2 primer start'
    elif primer_ends_before_seq_start(alignment, hxb2_segment_end):
        result['error'] = f'{end} sample probe starts after hxb2 primer end'

    # Get the start coordinate of the primer within the sample segment relative to hxb2
    primer_in_sample_probe_hxb2_start_coord = max(alignment.start,
                                                  primer['hxb2_start'])
    primer_in_sample_probe_hxb2_end_coord = min(alignment.start + probelen,
                                                primer['hxb2_end'])


def test_mismatch_revprimer(data):
    seq = None
    name = 'rev'
    primerbuffer = 20
    # probelen = 100
    # with open(data['cwd'] / 'data' / 'mismatched_revprimer.txt') as f:
    with open(data['cwd'] / 'data' / 'real_example.txt') as f:
        seq = f.read().strip()
    hxb2_target_start = primers[name]['hxb2_start'] - primerbuffer
    hxb2_target_end = primers[name]['hxb2_end']
    hxb2_target_seq = hxb2[hxb2_target_start:hxb2_target_end]
    print(hxb2_target_seq)
    print(data)
    # truncate seq to appropriate size
    # seq = seq[-probelen:]

    # finder = ProbeFinder(hxb2_target_seq, seq)
    finder = ProbeFinder(seq, hxb2_target_seq)
    print(finder)
    finder.start += hxb2_target_start
    print(finder)
    validated = validate_primer(finder, seq, primers[name])
    print(validated)