"""
Test suite for the seqlen column formatting fix in outcome_summary.py.

This module tests that integer seqlen values are written as integers,
not floats (e.g., "5" instead of "5.0").
"""
import csv
import io
import os
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from cfeproviral.outcome_summary import OutcomeSummary


class TestSeqlenFormatting:
    """Test that seqlen values are written as integers, not floats."""

    def test_seqlen_integer_formatting(self):
        """Test that integer seqlen values are written without decimal points."""
        # Create mock data with float seqlen values (as would come from pandas)
        conseqs_data = {
            'sample': ['a_S1', 'b_S2'],
            'run_name': ['test_conseqs', 'test_conseqs'],
            'seqtype': ['conseqs', 'conseqs'],
            'seqlen': [9016.0, 8885.0],  # Float values that should become integers
            'reference': ['HIV1-B-FR-K03455-seed', 'HIV1-B-FR-K03455-seed'],
            'sequence': ['ATCG' * 2254, 'ATCG' * 2221],
            'is_rev_comp': ['N', 'N'],
            'error': ['', ''],
            'fwd_error': ['', ''],
            'rev_error': ['', '']
        }
        
        contigs_data = {
            'sample': ['c_S3'],
            'run_name': ['test_contigs'],
            'seqtype': ['contigs'],
            'seqlen': [4419.0],  # Float value that should become integer
            'reference': ['HIV1-B-FR-K03455-seed'],
            'sequence': ['ATCG' * 1105],
            'is_rev_comp': ['N'],
            'error': [''],
            'fwd_error': [''],
            'rev_error': ['']
        }
        
        conseqs_df = pd.DataFrame(conseqs_data)
        contigs_df = pd.DataFrame(contigs_data)
        
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            outpath = Path(temp_dir)
            
            # Create the OutcomeSummary which should process and write the data
            summary = OutcomeSummary(conseqs_df, contigs_df, outpath, force_all_proviral=True)
            
            # Read the generated CSV file
            outcome_file = outpath / 'outcome_summary.csv'
            assert outcome_file.exists()
            
            with open(outcome_file, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                rows = list(reader)
            
            # Check that seqlen values are written as integers, not floats
            assert len(rows) == 3
            
            # Find the rows by sample name
            a_s1_row = next(row for row in rows if row['sample'] == 'a_S1')
            b_s2_row = next(row for row in rows if row['sample'] == 'b_S2')
            c_s3_row = next(row for row in rows if row['sample'] == 'c_S3')
            
            # Assert that seqlen values don't have decimal points
            assert a_s1_row['seqlen'] == '9016', f"Expected '9016', got '{a_s1_row['seqlen']}'"
            assert b_s2_row['seqlen'] == '8885', f"Expected '8885', got '{b_s2_row['seqlen']}'"
            assert c_s3_row['seqlen'] == '4419', f"Expected '4419', got '{c_s3_row['seqlen']}'"

    def test_seqlen_with_empty_values(self):
        """Test that empty seqlen values are handled correctly."""
        # Create mock data with empty seqlen values
        conseqs_data = {
            'sample': ['a_S1'],
            'run_name': ['test_conseqs'],
            'seqtype': ['conseqs'],
            'seqlen': [''],  # Empty string
            'reference': [''],
            'sequence': [''],
            'is_rev_comp': [''],
            'error': ['no sequence'],
            'fwd_error': [''],
            'rev_error': ['']
        }
        
        contigs_data = {
            'sample': [],
            'run_name': [],
            'seqtype': [],
            'seqlen': [],
            'reference': [],
            'sequence': [],
            'is_rev_comp': [],
            'error': [],
            'fwd_error': [],
            'rev_error': []
        }
        
        conseqs_df = pd.DataFrame(conseqs_data)
        contigs_df = pd.DataFrame(contigs_data)
        
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            outpath = Path(temp_dir)
            
            # Create the OutcomeSummary which should process and write the data
            summary = OutcomeSummary(conseqs_df, contigs_df, outpath, force_all_proviral=True)
            
            # Read the generated CSV file
            outcome_file = outpath / 'outcome_summary.csv'
            assert outcome_file.exists()
            
            with open(outcome_file, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                rows = list(reader)
            
            # Check that empty seqlen values remain empty
            assert len(rows) == 1
            assert rows[0]['seqlen'] == ''

    def test_seqlen_with_nan_values(self):
        """Test that NaN seqlen values are handled correctly."""
        # Create mock data with NaN seqlen values
        conseqs_data = {
            'sample': ['a_S1'],
            'run_name': ['test_conseqs'],
            'seqtype': ['conseqs'],
            'seqlen': [float('nan')],  # NaN value
            'reference': [''],
            'sequence': [''],
            'is_rev_comp': [''],
            'error': ['no sequence'],
            'fwd_error': [''],
            'rev_error': ['']
        }
        
        contigs_data = {
            'sample': [],
            'run_name': [],
            'seqtype': [],
            'seqlen': [],
            'reference': [],
            'sequence': [],
            'is_rev_comp': [],
            'error': [],
            'fwd_error': [],
            'rev_error': []
        }
        
        conseqs_df = pd.DataFrame(conseqs_data)
        contigs_df = pd.DataFrame(contigs_data)
        
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            outpath = Path(temp_dir)
            
            # Create the OutcomeSummary which should process and write the data
            summary = OutcomeSummary(conseqs_df, contigs_df, outpath, force_all_proviral=True)
            
            # Read the generated CSV file
            outcome_file = outpath / 'outcome_summary.csv'
            assert outcome_file.exists()
            
            with open(outcome_file, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                rows = list(reader)
            
            # Check that NaN seqlen values are preserved as-is (should not cause errors)
            assert len(rows) == 1
            # NaN values get converted to string representation, so check for that
            assert 'nan' in rows[0]['seqlen'].lower() or rows[0]['seqlen'] == ''

    def test_seqlen_with_string_numeric_values(self):
        """Test that string numeric seqlen values are converted to integers."""
        # Create mock data with string seqlen values
        conseqs_data = {
            'sample': ['a_S1', 'b_S2'],
            'run_name': ['test_conseqs', 'test_conseqs'],
            'seqtype': ['conseqs', 'conseqs'],
            'seqlen': ['9016.0', '8885'],  # String values that should become integers
            'reference': ['HIV1-B-FR-K03455-seed', 'HIV1-B-FR-K03455-seed'],
            'sequence': ['ATCG' * 2254, 'ATCG' * 2221],
            'is_rev_comp': ['N', 'N'],
            'error': ['', ''],
            'fwd_error': ['', ''],
            'rev_error': ['', '']
        }
        
        contigs_data = {
            'sample': [],
            'run_name': [],
            'seqtype': [],
            'seqlen': [],
            'reference': [],
            'sequence': [],
            'is_rev_comp': [],
            'error': [],
            'fwd_error': [],
            'rev_error': []
        }
        
        conseqs_df = pd.DataFrame(conseqs_data)
        contigs_df = pd.DataFrame(contigs_data)
        
        # Create a temporary directory for output
        with tempfile.TemporaryDirectory() as temp_dir:
            outpath = Path(temp_dir)
            
            # Create the OutcomeSummary which should process and write the data
            summary = OutcomeSummary(conseqs_df, contigs_df, outpath, force_all_proviral=True)
            
            # Read the generated CSV file
            outcome_file = outpath / 'outcome_summary.csv'
            assert outcome_file.exists()
            
            with open(outcome_file, 'r', newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                rows = list(reader)
            
            # Check that seqlen values are written as integers
            assert len(rows) == 2
            
            sample1_row = next(row for row in rows if row['sample'] == 'a_S1')
            sample2_row = next(row for row in rows if row['sample'] == 'b_S2')
            
            assert sample1_row['seqlen'] == '9016', f"Expected '9016', got '{sample1_row['seqlen']}'"
            assert sample2_row['seqlen'] == '8885', f"Expected '8885', got '{sample2_row['seqlen']}'"
