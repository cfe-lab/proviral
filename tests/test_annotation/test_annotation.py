import os
from pathlib import Path

import cfeproviral.utils as utils


def test_annotation():
    cwd = Path(os.path.realpath(__file__)).parent
    actual_annot = utils.mod_annot
    expected_annot = {x['gene']: [int(x['start']), int(x['stop'])]
                      for x in utils.read_csv(cwd / 'valid' / 'valid_mod_annot.csv')}
    utils.csv_to_bed(cwd / 'valid' / 'valid_mod_annot.csv', 'MOD_HXB2')
    bed_path: Path = cwd / 'valid' / 'valid_mod_annot.csv.bed'
    bed_path.unlink()

    assert actual_annot == expected_annot
