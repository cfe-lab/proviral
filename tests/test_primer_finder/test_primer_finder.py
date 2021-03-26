import pytest
import math
import os
import sys
from pathlib import Path

cwd = Path(os.path.realpath(__file__)).parent

import gene_splicer.utils as utils
from gene_splicer.primer_finder_class import PrimerFinder
import Levenshtein
from gene_splicer.primer_finder import primers


@pytest.fixture
def data():
    return {'cwd': Path(os.path.realpath(__file__)).parent}


def test_case():
    case = {'sequence': 'AAA'}

    sample_size = 50
    extended_length = 200

    fwd_primer = PrimerFinder(case['sequence'],
                              primers['fwd']['seq'],
                              'fwd',
                              primers['fwd']['hxb2_start'],
                              primers['fwd']['hxb2_end'],
                              sample_size=sample_size)

    if not fwd_primer.is_full_length:
        fwd_primer2 = PrimerFinder(case['sequence'],
                                   primers['fwd']['seq'],
                                   'fwd',
                                   primers['fwd']['hxb2_start'],
                                   primers['fwd']['hxb2_end'],
                                   sample_size=extended_length)
        if fwd_primer2.is_full_length:
            fwd_primer = fwd_primer2

    rev_primer = PrimerFinder(case['sequence'],
                              primers['rev']['seq'],
                              'rev',
                              primers['rev']['hxb2_start'],
                              primers['rev']['hxb2_end'],
                              sample_size=sample_size)

    if not rev_primer.is_full_length:
        rev_primer2 = PrimerFinder(case['sequence'],
                                   primers['rev']['seq'],
                                   'rev',
                                   primers['rev']['hxb2_start'],
                                   primers['rev']['hxb2_end'],
                                   sample_size=extended_length)
        if rev_primer2.is_full_length:
            rev_primer = rev_primer2

    print()
    print(fwd_primer)
    print(rev_primer)


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