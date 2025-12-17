"""
Tests for CFEIntact output version column post-processing.
"""
import csv
from cfeproviral.primer_finder import add_versions_to_cfeintact_outputs


def test_add_versions_to_cfeintact_outputs_with_qseqid(tmp_path):
    """Test that version columns are added correctly when qseqid is present."""
    # Create a dummy CFEIntact CSV file with qseqid
    csv_file = tmp_path / "holistic.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['qseqid', 'verdict', 'score'])
        writer.writeheader()
        writer.writerow({'qseqid': 'seq1::sample1::HXB2::conseq', 'verdict': 'Intact', 'score': '100'})
        writer.writerow({'qseqid': 'seq2::sample2::HXB2::contig', 'verdict': 'Defect', 'score': '50'})

    # Mock all_samples dictionary
    all_samples = {
        'sample1': {'micall_version': '7.15.0'},
        'sample2': {'micall_version': '7.16.1'}
    }

    # Call the function
    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    # Verify the output
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # Check if new columns are present
    assert 'cfeproviral_version' in fieldnames
    assert 'cfeintact_version' in fieldnames
    assert 'micall_version' in fieldnames

    # Check values
    assert rows[0]['micall_version'] == '7.15.0'
    assert rows[1]['micall_version'] == '7.16.1'
    
    # Check that original columns are preserved
    assert rows[0]['verdict'] == 'Intact'
    assert rows[0]['score'] == '100'
    assert rows[1]['verdict'] == 'Defect'
    assert rows[1]['score'] == '50'

    # Check that versions are not empty
    assert rows[0]['cfeproviral_version']
    assert rows[0]['cfeintact_version']
    assert rows[1]['cfeproviral_version']
    assert rows[1]['cfeintact_version']


def test_add_versions_to_cfeintact_outputs_without_qseqid(tmp_path):
    """Test that version columns are added even without qseqid, with None for micall_version."""
    # Create a CSV file without qseqid column
    csv_file = tmp_path / "summary.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['metric', 'value'])
        writer.writeheader()
        writer.writerow({'metric': 'total_sequences', 'value': '100'})
        writer.writerow({'metric': 'intact_count', 'value': '75'})

    # Mock all_samples dictionary (won't be used for this file)
    all_samples = {
        'sample1': {'micall_version': '7.15.0'}
    }

    # Call the function
    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    # Verify the output
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # Check if new columns are present
    assert 'cfeproviral_version' in fieldnames
    assert 'cfeintact_version' in fieldnames
    assert 'micall_version' in fieldnames

    # Check that micall_version is empty/None when qseqid is not present
    assert not rows[0]['micall_version'] or rows[0]['micall_version'] == 'None'
    assert not rows[1]['micall_version'] or rows[1]['micall_version'] == 'None'

    # Check that other versions are present
    assert rows[0]['cfeproviral_version']
    assert rows[0]['cfeintact_version']

    # Check that original data is preserved
    assert rows[0]['metric'] == 'total_sequences'
    assert rows[0]['value'] == '100'


def test_add_versions_to_cfeintact_outputs_malformed_qseqid(tmp_path):
    """Test handling of malformed qseqid that doesn't follow name::sample::ref::type format."""
    csv_file = tmp_path / "defects.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['qseqid', 'defect_type'])
        writer.writeheader()
        writer.writerow({'qseqid': 'invalid_format', 'defect_type': 'APOBEC'})
        writer.writerow({'qseqid': 'only::one::colon', 'defect_type': 'PSI'})
        writer.writerow({'qseqid': 'name::sample::ref::type', 'defect_type': 'Stop'})

    all_samples = {
        'sample': {'micall_version': '7.15.0'}
    }

    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Malformed qseqid should result in None for micall_version
    assert not rows[0]['micall_version'] or rows[0]['micall_version'] == 'None'
    assert not rows[1]['micall_version'] or rows[1]['micall_version'] == 'None'
    # Valid qseqid should retrieve the version
    assert rows[2]['micall_version'] == '7.15.0'

    # But other versions should always be present
    for row in rows:
        assert row['cfeproviral_version']
        assert row['cfeintact_version']


def test_add_versions_to_cfeintact_outputs_missing_sample(tmp_path):
    """Test handling when sample is not in all_samples dictionary."""
    csv_file = tmp_path / "output.csv"
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['qseqid', 'data'])
        writer.writeheader()
        writer.writerow({'qseqid': 'name::unknown_sample::ref::type', 'data': 'test'})

    all_samples = {
        'known_sample': {'micall_version': '7.15.0'}
    }

    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Sample not in dictionary should result in None
    assert not rows[0]['micall_version'] or rows[0]['micall_version'] == 'None'
    
    # But other versions should be present
    assert rows[0]['cfeproviral_version']
    assert rows[0]['cfeintact_version']


def test_add_versions_to_cfeintact_outputs_multiple_files(tmp_path):
    """Test that all CSV files in the directory are processed."""
    # Create multiple CSV files
    for filename in ['holistic.csv', 'defects.csv', 'genes.csv']:
        csv_file = tmp_path / filename
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['qseqid', 'col1'])
            writer.writeheader()
            writer.writerow({'qseqid': 'name::sample1::ref::type', 'col1': 'value'})

    all_samples = {
        'sample1': {'micall_version': '7.15.0'}
    }

    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    # Verify all files were processed
    for filename in ['holistic.csv', 'defects.csv', 'genes.csv']:
        csv_file = tmp_path / filename
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            rows = list(reader)

        assert 'cfeproviral_version' in fieldnames
        assert 'cfeintact_version' in fieldnames
        assert 'micall_version' in fieldnames
        assert rows[0]['micall_version'] == '7.15.0'


def test_add_versions_to_cfeintact_outputs_empty_file(tmp_path):
    """Test handling of empty CSV files."""
    csv_file = tmp_path / "empty.csv"
    with open(csv_file, 'w', newline='') as f:
        f.write('')

    all_samples = {}

    # Should not crash on empty file
    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    # File should remain empty or minimally modified
    assert csv_file.exists()


def test_add_versions_to_cfeintact_outputs_preserves_column_order(tmp_path):
    """Test that version columns are added at the end, preserving original order."""
    csv_file = tmp_path / "test.csv"
    original_fieldnames = ['qseqid', 'col_a', 'col_b', 'col_c']
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=original_fieldnames)
        writer.writeheader()
        writer.writerow({'qseqid': 'name::sample1::ref::type', 'col_a': '1', 'col_b': '2', 'col_c': '3'})

    all_samples = {
        'sample1': {'micall_version': '7.15.0'}
    }

    add_versions_to_cfeintact_outputs(tmp_path, all_samples)

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames)

    # Check that original columns come first
    assert fieldnames[:4] == original_fieldnames
    # Check that version columns are added at the end
    assert fieldnames[4:] == ['cfeproviral_version', 'cfeintact_version', 'micall_version']
