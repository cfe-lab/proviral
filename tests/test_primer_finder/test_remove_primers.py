
"""
Verify that remove_primers() correctly chops off only the synthetic
forward- and reverse-primer "handles" and returns exactly the payload.

Each test constructs a fake contig by embedding forward and/or reverse
primer fragments around a payload, runs the detection and trimming steps
(`find_primers` → `filter_df` → `remove_primers`), and asserts that trimming
behavior matches expectations (or that no trimming occurs when primers
overlap or are missing).

Caveat: the pipeline’s QC filter will drop any sequence containing
non-TCGA characters (e.g. digits); such rows are silently removed before
trimming, which can make tests “disappear.” If your payload includes
non-TCGA bases, either clear its resulting error fields or restrict payloads
to A/C/T/G only.

Maintenance: if primer definitions or QC rules change, these tests
will need to be updated accordingly.
"""

import pandas as pd
import numpy as np
import pytest
import re
from cfeproviral.primer_finder_class import PrimerFinder
from cfeproviral.primer_finder import remove_primers, filter_df, primers, PrimerFinderErrors, record_primers, find_primers


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


def test_remove_primers_is_applied_by_filter_df():
    # --- construct a fake molecule -------------------------------
    payload       = "PAYLOAD"
    fwd_size      = 3
    fwd_start     = 2            # forward primer begins after 2 bases of pad
    rev_size      = 4            # reverse primer length
    front_pad     = "X" * fwd_start
    fwd_primer    = "F" * fwd_size
    rev_primer    = "R" * rev_size
    back_pad      = "Y" * 5      # some extra tail padding

    # we must compute absolute index where the reverse primer starts:
    # front_pad + fwd_primer + payload
    rev_start = fwd_start + fwd_size + len(payload)

    # assemble the full sequence:
    full_seq = front_pad + fwd_primer + payload + rev_primer + back_pad

    # sanity:
    assert full_seq[ rev_start : rev_start + rev_size ] == rev_primer

    # --- build a DataFrame row just like primer_finder leaves it ----
    sample_size = len(full_seq)
    df = pd.DataFrame([{
        # must be “no error” so filter_df doesn’t drop the row
        "error":        np.nan,
        "fwd_error":    np.nan,
        "rev_error":    np.nan,
        # any reference that doesn’t contain “reverse” or “unknown”
        "reference":    "sample1",
        "seqtype":      "contig",
        # the long sequence, to be trimmed
        "sequence":     full_seq,
        # the four metadata fields remove_primers needs:
        "fwd_sample_primer_start": fwd_start,
        "fwd_sample_primer_size":  fwd_size,
        "rev_sample_primer_start": rev_start,
        "rev_sample_primer_size":  rev_size,
    }])

    # run filter_df without deduplication
    filtered = filter_df(sample_size, df, nodups=False)

    # after trimming, we should have exactly one row…
    assert len(filtered) == 1

    # …and its “sequence” column must be exactly the payload
    assert filtered.iloc[0]["sequence"] == payload


@pytest.mark.parametrize("direction", ("fwd","rev"))
def test_primerfinder_detects_pure_acgt_suffixes(direction):
    # grab the canonical primer (may contain mixtures Y/R/etc)
    primer_seq = primers[direction]["seq"]
    other      = primers["rev" if direction=="fwd" else "fwd"]["seq"]

    # embed it in a little fake molecule
    if direction == "fwd":
        full = primer_seq + "PAYLOAD" + other
    else:
        full = other + "PAYLOAD" + primer_seq

    sample_size     = len(full)
    validation_size = 5

    finder = PrimerFinder(
        full_sample     = full,
        primer          = primer_seq,
        direction       = direction,
        hxb2_start      = primers[direction]["hxb2_start"],
        hxb2_end        = primers[direction]["hxb2_end"],
        validation_size = validation_size,
        sample_size     = sample_size,
    )

    # 1) primer_primer is pure A/C/G/T
    found = finder.sample_primer
    assert found and all(nt in "ACGT" for nt in found)

    # 2) it must match the longest A/C/G/T suffix of the original
    pure_tail = re.search(r"([ACGT]+)$", primer_seq).group(1)
    assert found == pure_tail

    # 3) only the reverse primer (no mixtures) is long enough to auto-validate
    if direction == "rev":
        assert finder.is_valid,        "reverse primer should be valid"
        assert finder.is_full_length,  "reverse primer should be full-length"
    else:
        assert not finder.is_valid,       "forward primer suffix <7bp so not valid"
        assert not finder.is_full_length, "forward suffix never reaches full_length"


def setup_trim_case(tmp_path, fwd_frag, rev_frag, payload,
                    front_pad='', back_pad='', sample_name="TEST_SAMPLE"):
    """Build a fake contig with given primer fragments and payload,
    run the primer‐finder pipeline, and return the trimmed DataFrame."""
    # 1) Make output dir
    outdir = tmp_path / "out"
    outdir.mkdir()

    # 2) Assemble full sequence
    full_seq = f"{front_pad}{fwd_frag}{payload}{rev_frag}{back_pad}"

    # 3) Write contigs.csv
    contigs_csv = tmp_path / "contigs.csv"
    contigs_csv.write_text(f"sample,region,sequence\n{sample_name},HIV_dummy,{full_seq}\n")

    # 4) Write cascade.csv so it's treated as proviral
    cascade_csv = tmp_path / "cascade.csv"
    cascade_csv.write_text(f"sample,remap\n{sample_name},1\n")

    # 5) Run primer-finder
    analysis_csv = find_primers(
        csv_filepath=contigs_csv.open("r"),
        outpath=outdir,
        run_name="run1",
        all_samples={sample_name: 1},
        seqtype="contigs",
        sample_size=len(full_seq),
        force_all_proviral=True
    )

    # 6) Load results and clear all error flags
    df = pd.read_csv(analysis_csv)
    for col in ("error", "fwd_error", "rev_error"):
        df[col] = pd.NA

    # 7) Add the seqtype column expected by filter_df
    df["seqtype"] = "contigs"

    # 8) Apply filter_df (which calls remove_primers under the hood)
    trimmed = filter_df(sample_size=len(full_seq), df=df, nodups=False)
    return trimmed, full_seq


def test_partial_primer_trimming(tmp_path):
    # Use partial‐primer fragments:
    fwd_frag = primers["fwd"]["nomix"][-10:]
    rev_frag = primers["rev"]["nomix"][:8]
    front_pad = "AAA"
    back_pad = "TTT"

    # Payload must be A/C/T/G only to avoid non-TCGA QC filtering
    payload = "ACTAAACTATATATTTAAATATATGTTCTCTATTC"

    trimmed, full_seq = setup_trim_case(
        tmp_path, fwd_frag, rev_frag, payload, front_pad, back_pad
    )

    # Should have exactly one row, and its sequence equals the payload
    assert len(trimmed) == 1
    assert trimmed.iloc[0]["sequence"] == payload
