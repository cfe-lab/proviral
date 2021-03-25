from io import StringIO
from pathlib import Path

from gene_splicer.study_summary import StudySummary


def test_load_runs():
    summary = StudySummary()
    requested_runs = ['/fake/path/200101_M11111_0001_0000000-J1HRQ',
                      '/fake/path/200108_M11111_0002_0000000-P3RS4']

    summary.load_runs(requested_runs)

    assert summary.run_paths == tuple(Path(run_path)
                                      for run_path in requested_runs)


def test_load_samples(tmp_path):
    runs_root = tmp_path
    run1_path = runs_root / '200101_M11111_0001_0000000-J1HRQ'
    run1_path.mkdir()
    run2_path = runs_root / '200115_M22222_0003_0000000-Y8E4T'
    run2_path.mkdir()

    samples_csv = StringIO("""\
sample,run,pid
E0001_S1,200101_M11111,P1
E0002_S2,200101_M11111,P2
E0003_S1,200115_M22222,P2
""")
    summary = StudySummary()

    summary.load_samples(samples_csv, runs_root)

    assert summary.run_paths == (run1_path, run2_path)


def test_load_samples_then_runs(tmp_path):
    runs_root = tmp_path
    run1_path = runs_root / '200101_M11111_0001_0000000-J1HRQ'
    run1_path.mkdir()
    run2_path = runs_root / '200115_M22222_0003_0000000-Y8E4T'
    run2_path.mkdir()

    samples_csv = StringIO("""\
sample,run,pid
E0001_S1,200101_M11111,P1
E0002_S2,200101_M11111,P2
E0003_S1,200115_M22222,P2
""")
    requested_runs = [str(run1_path)]
    summary = StudySummary()

    summary.load_samples(samples_csv, runs_root)
    summary.load_runs(requested_runs)

    assert summary.run_paths == (run1_path,)


def test_load_runs_then_samples(tmp_path):
    runs_root = tmp_path
    run1_path = runs_root / '200101_M11111_0001_0000000-J1HRQ'
    run1_path.mkdir()
    run2_path = runs_root / '200115_M22222_0003_0000000-Y8E4T'
    run2_path.mkdir()

    samples_csv = StringIO("""\
sample,run,pid
E0001-NFLHIVDNA_S1,200101_M11111,P1
E0002-NFLHIVDNA_S2,200101_M11111,P2
E0003-NFLHIVDNA_S1,200115_M22222,P2
""")
    requested_runs = [str(run1_path)]
    summary = StudySummary()

    summary.load_runs(requested_runs)
    summary.load_samples(samples_csv, runs_root)

    assert summary.run_paths == (run1_path,)


def test_load_outcome():
    summary = StudySummary()
    outcome_summary_csv = StringIO("""\
sample,run,passed,error
P1-NFLHIVDNA_S1,200101_M11111,True,
P2-NFLHIVDNA_S2,200101_M11111,False,primer error
P3-NFLHIVDNA_S3,200101_M11111,False,multiple contigs
""")
    expected_run_counts = {'200101_M11111': dict(samples=3,
                                                 passed=1,
                                                 errors=2,
                                                 primer_error=1,
                                                 multiple_contigs=1)}

    summary.load_outcome(outcome_summary_csv)

    assert summary.run_counts == expected_run_counts


def test_load_multiple_run_outcomes():
    summary = StudySummary()
    outcome1_summary_csv = StringIO("""\
sample,run,passed,error
P1-NFLHIVDNA_S1,200101_M11111,True,
P2-NFLHIVDNA_S2,200101_M11111,False,primer error
""")
    outcome2_summary_csv = StringIO("""\
sample,run,passed,error
P2-NFLHIVDNA_S1,200108_M11111,True,
P3-NFLHIVDNA_S2,200108_M11111,False,multiple contigs
""")
    expected_run_counts = {'200101_M11111': dict(samples=2,
                                                 passed=1,
                                                 errors=1,
                                                 primer_error=1),
                           '200108_M11111': dict(samples=2,
                                                 passed=1,
                                                 errors=1,
                                                 multiple_contigs=1)}
    expected_participant_counts = {'P1': dict(samples=1, passed=1),
                                   'P2': dict(samples=2,
                                              passed=1,
                                              errors=1,
                                              primer_error=1),
                                   'P3': dict(samples=1,
                                              errors=1,
                                              multiple_contigs=1)}

    summary.load_outcome(outcome1_summary_csv)
    summary.load_outcome(outcome2_summary_csv)

    assert summary.run_counts == expected_run_counts
    assert summary.participant_counts == expected_participant_counts


def test_participant_lookup(tmp_path):
    runs_root = tmp_path
    run1_path = runs_root / '200101_M11111_0001_0000000-J1HRQ'
    run1_path.mkdir()

    samples_csv = StringIO("""\
sample,run,pid
E0001-NFLHIVDNA_S1,200101_M11111,P1
E0002-NFLHIVDNA_S2,200101_M11111,P2
E0003-NFLHIVDNA_S3,200101_M11111,P2
""")
    outcome_summary_csv = StringIO("""\
sample,run,passed,error
E0001-NFLHIVDNA_S1,200101_M11111,True,
E0002-NFLHIVDNA_S2,200101_M11111,False,primer error
E0003-NFLHIVDNA_S3,200101_M11111,False,multiple contigs
""")
    expected_participant_counts = {'P1': dict(samples=1, passed=1),
                                   'P2': dict(samples=2,
                                              errors=2,
                                              primer_error=1,
                                              multiple_contigs=1)}
    summary = StudySummary()

    summary.load_samples(samples_csv, runs_root)
    summary.load_outcome(outcome_summary_csv)

    assert summary.participant_counts == expected_participant_counts


def test_write():
    summary = StudySummary()
    outcome_summary_csv = StringIO("""\
sample,run,passed,error
P1-NFLHIVDNA_S1,200101_M11111,True,
P2-NFLHIVDNA_S2,200101_M11111,False,primer error
P3-NFLHIVDNA_S3,200101_M11111,False,multiple contigs
""")
    expected_study_summary_csv = """\
type,name,samples,passed,errors,no_sequence,non_hiv,hiv_but_failed,\
primer_error,low_internal_cov,multiple_contigs
run,200101_M11111,3,1,2,0,0,0,1,0,1
participant,P1,1,1,0,0,0,0,0,0,0
participant,P2,1,0,1,0,0,0,1,0,0
participant,P3,1,0,1,0,0,0,0,0,1
total,total,3,1,2,0,0,0,1,0,1
"""
    study_summary_csv = StringIO()

    summary.load_outcome(outcome_summary_csv)
    summary.write(study_summary_csv)

    assert study_summary_csv.getvalue() == expected_study_summary_csv
