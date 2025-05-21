
import pytest
from cfeproviral.primer_finder import remove_primers

class DummyRow:
    """
    A minimal stand-in with exactly the 5 attributes remove_primers uses:
      .sequence                   (the full string)
      .fwd_sample_primer_start    (int)
      .fwd_sample_primer_size     (int)
      .rev_sample_primer_start    (int)
      .rev_sample_primer_size     (int)
    """
    def __init__(self,
                 full_sequence: str,
                 fwd_start: int,
                 fwd_size: int,
                 rev_start: int,
                 rev_size: int):
        self.sequence = full_sequence
        self.fwd_sample_primer_start = fwd_start
        self.fwd_sample_primer_size  = fwd_size
        self.rev_sample_primer_start = rev_start
        self.rev_sample_primer_size  = rev_size

@pytest.mark.parametrize(
    "payload, fwd_start, fwd_size, rev_size",
    [
        # primers flush exactly to the edges
        ("PAYLOAD",      0, 3, 3),
        # forward primer offset 2 bases in, reverse flush
        ("HELLOWORLD",   2, 4, 3),
        # forward flush, reverse primer length 5
        ("ABCDEFG",      0, 5, 4),
        # both ends offset
        ("1234567890",   1, 3, 2),
        # zero-length forward primer
        ("X",            0, 0, 1),
        # zero-length reverse primer
        ("Z",            1, 1, 0),
    ]
)
def test_remove_primers_recovers_payload(payload, fwd_start, fwd_size, rev_size):
    # 1) build the forward primer and payload and reverse primer:
    fwd_primer = "F" * fwd_size
    rev_primer = "R" * rev_size

    # 2) pad the front so that fwd_sample_primer_start==fwd_start
    front_pad = "X" * fwd_start

    # 3) assemble full molecule
    full = front_pad + fwd_primer + payload + rev_primer

    # 4) compute the absolute index of the reverse primer
    #    (i.e. where 'R'*rev_size starts in full)
    rev_abs = len(front_pad) + fwd_size + len(payload)

    # 5) sanity check
    assert full[rev_abs: rev_abs + rev_size] == rev_primer

    # 6) call remove_primers with sample_size == len(full)
    row = DummyRow(full,
                   fwd_start=fwd_start,
                   fwd_size=fwd_size,
                   rev_start=rev_abs,
                   rev_size=rev_size)

    out = remove_primers(len(full), row)

    # should return the same object, and strip off exactly the primers:
    assert out is row
    assert row.sequence == payload


def test_remove_primers_no_trimming_if_all_zero():
    """
    If both primer-sizes are zero, remove_primers should not touch .sequence.
    """
    payload = "SOME_SEQ"
    # full == payload, no primers at all
    row = DummyRow(payload, fwd_start=0, fwd_size=0, rev_start=0, rev_size=0)
    out = remove_primers(len(payload), row)
    assert out.sequence == payload
