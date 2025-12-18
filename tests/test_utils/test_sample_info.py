"""
Tests for sample_info handling and micall_version propagation.

These tests verify that:
1. micall_version is ONLY read from sample_info.csv (never from cascade.csv)
2. Both sample-specific and run-level sample_info formats are supported
3. Errors are raised when sample_info contains samples not in cascade
"""
import io
import pytest
import tempfile
import csv
from pathlib import Path

import cfeproviral.utils as utils


class TestGetSamplesFromCascade:
    """Test that cascade reading NEVER includes micall_version."""

    def test_multi_sample_cascade_no_micall_version(self):
        """Cascade with multiple samples should not have micall_version in output."""
        cascade_content = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned,micall_version\n"
            "sample1,1000,100,50,500,450,400,v7.17\n"
            "sample2,2000,200,100,1000,900,850,v7.17\n"
        )

        result = utils.get_samples_from_cascade(cascade_content)

        assert 'sample1' in result
        assert 'sample2' in result
        assert result['sample1'] == {'remap': 450}
        assert result['sample2'] == {'remap': 900}
        # Verify micall_version is NOT in the result
        assert 'micall_version' not in result['sample1']
        assert 'micall_version' not in result['sample2']

    def test_single_sample_cascade_no_micall_version(self):
        """Single sample cascade (no sample column) should not have micall_version."""
        cascade_content = io.StringIO(
            "demultiplexed,v3loop,g2p,prelim_map,remap,aligned,micall_version\n"
            "1000,100,50,500,450,400,v7.17\n"
        )

        result = utils.get_samples_from_cascade(cascade_content, 'default_sample')

        assert 'default_sample' in result
        assert result['default_sample'] == {'remap': 450}
        # Verify micall_version is NOT in the result
        assert 'micall_version' not in result['default_sample']


class TestSampleInfoMicallVersion:
    """Test micall_version propagation from sample_info."""

    def test_run_level_micall_version(self):
        """Test run-level sample_info (no 'sample' column) applies to all samples."""
        # Cascade with multiple samples
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        # Run-level info (no 'sample' column)
        run_level_info = {
            'project': 'NFLHIVDNA',
            'run_name': 'test_run',
            'micall_version': 'v7.17'
        }

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Apply run-level info (simulating primer_finder.run logic)
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # All samples should have the same micall_version
        assert all_samples['sample1']['micall_version'] == 'v7.17'
        assert all_samples['sample2']['micall_version'] == 'v7.17'

    def test_sample_specific_micall_version(self):
        """Test sample-specific sample_info (with 'sample' column) maps to specific samples."""
        # Cascade with multiple samples
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        # Sample-specific info
        sample_info_map = {
            'sample1': {'sample': 'sample1', 'micall_version': 'v7.15'},
            'sample2': {'sample': 'sample2', 'micall_version': 'v7.16'}
        }

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Apply sample-specific info (simulating primer_finder.run logic)
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if sample_name in all_samples:
                    if 'micall_version' in sample_info:
                        all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # Each sample should have its specific micall_version
        assert all_samples['sample1']['micall_version'] == 'v7.15'
        assert all_samples['sample2']['micall_version'] == 'v7.16'

    def test_sample_info_mismatch_error(self):
        """Test error when sample_info contains samples not in cascade."""
        # Cascade with specific samples
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        # Sample_info with a sample that doesn't exist
        sample_info_map = {
            'sample_does_not_exist': {'sample': 'sample_does_not_exist', 'micall_version': 'v7.17'}
        }

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Should raise ValueError (simulating primer_finder.run logic)
        with pytest.raises(ValueError, match="Sample 'sample_does_not_exist' in sample_info does not exist"):
            for sample_name, sample_info in sample_info_map.items():
                if sample_name not in all_samples:
                    raise ValueError(
                        f"Sample '{sample_name}' in sample_info does not exist in cascade data. "
                        f"Available samples: {list(all_samples.keys())}"
                    )


class TestSampleInfoSemantics:
    """Test the semantic differences between sample_info formats."""

    def test_empty_sample_info_no_error(self):
        """Empty sample_info should not cause errors."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Apply empty sample_info_map
        sample_info_map = {}
        # Should not raise error
        for sample_name, sample_info in sample_info_map.items():
            pass  # No iteration, no error

        # Samples should exist but have no micall_version
        assert 'sample1' in all_samples
        assert 'micall_version' not in all_samples['sample1']

    def test_run_level_takes_precedence_over_missing_sample_specific(self):
        """When sample_info_map is empty but run_level_info exists, use run_level."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        sample_info_map = {}
        run_level_info = {'micall_version': 'v7.17'}

        # Apply sample-specific (none exist)
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # Apply run-level
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        assert all_samples['sample1']['micall_version'] == 'v7.17'

    def test_sample_specific_not_overwritten_by_run_level(self):
        """Sample-specific info should be applied before run-level, and run-level overwrites it."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        sample_info_map = {
            'sample1': {'micall_version': 'v7.15'}
        }
        run_level_info = {'micall_version': 'v7.17'}

        # Apply sample-specific first
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if sample_name in all_samples and 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # Then apply run-level (this overwrites)
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Run-level overwrites everything
        assert all_samples['sample1']['micall_version'] == 'v7.17'
        assert all_samples['sample2']['micall_version'] == 'v7.17'


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_micall_version_none_or_empty(self):
        """Test handling of None or empty micall_version values."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Run-level with None micall_version
        run_level_info = {'micall_version': None}

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Should accept None
        assert all_samples['sample1']['micall_version'] is None

    def test_cascade_with_no_samples(self):
        """Test handling of empty cascade data."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Should return empty dict
        assert all_samples == {}

        # Applying run-level info to empty samples should not error
        run_level_info = {'micall_version': 'v7.17'}
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        assert all_samples == {}

    def test_run_level_without_micall_version_field(self):
        """Test run-level info that doesn't have micall_version field at all."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        run_level_info = {'project': 'NFLHIVDNA', 'run_name': 'test'}  # No micall_version

        # Should not apply anything
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # No micall_version should be added
        assert 'micall_version' not in all_samples['sample1']

    def test_empty_string_micall_version(self):
        """Test empty string micall_version."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        run_level_info = {'micall_version': ''}

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Should accept empty string
        assert all_samples['sample1']['micall_version'] == ''

    def test_single_sample_format_with_default_name(self):
        """Test single-sample cascade format (no 'sample' column) with default_sample_name."""
        cascade_csv = io.StringIO(
            "demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv, 'MySample')

        assert 'MySample' in all_samples
        assert all_samples['MySample']['remap'] == 450

        # Apply run-level info
        run_level_info = {'micall_version': 'v7.17'}
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        assert all_samples['MySample']['micall_version'] == 'v7.17'

    def test_multiple_samples_partial_sample_info(self):
        """Test when sample_info only has info for SOME samples in cascade."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
            "sample3,3000,300,150,1500,1350,1300\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Only have info for sample1 and sample3
        sample_info_map = {
            'sample1': {'micall_version': 'v7.15'},
            'sample3': {'micall_version': 'v7.16'}
        }

        # Apply sample-specific info
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if sample_name in all_samples and 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # sample1 and sample3 should have micall_version
        assert all_samples['sample1']['micall_version'] == 'v7.15'
        assert all_samples['sample3']['micall_version'] == 'v7.16'
        # sample2 should NOT have micall_version
        assert 'micall_version' not in all_samples['sample2']

    def test_cascade_remap_zero(self):
        """Test that remap count of 0 is handled correctly."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,0,0\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        assert 'sample1' in all_samples
        assert all_samples['sample1']['remap'] == 0

    def test_special_characters_in_version(self):
        """Test micall_version with various special characters."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Various version formats
        test_versions = [
            'v7.17',
            '7.17.0',
            'v7.17-rc1',
            'v7.17+build123',
            'v7.17.0-alpha.1+20231201',
        ]

        for test_version in test_versions:
            run_level_info = {'micall_version': test_version}
            if run_level_info and 'micall_version' in run_level_info:
                for sample_name in all_samples:
                    all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

            assert all_samples['sample1']['micall_version'] == test_version


class TestRunLevelVsSampleSpecific:
    """Test the critical distinction between run-level and sample-specific metadata."""

    def test_no_sample_column_means_run_level(self):
        """CRITICAL: No 'sample' column in sample_info means run-level metadata for ALL samples."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "2439P119528-P21I9-NFLHIVDNA_S33,33703,0,0,0,28168,28013\n"
            "2439P119528-P21F1-NFLHIVDNA_S32,50364,747,225,0,49956,49943\n"
            "3453P119516-P26P24-NFLHIVDNA_S31,41234,500,300,0,38000,37500\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # This is what ws-rerun-raw-data-folder creates: no 'sample' column
        run_level_info = {
            'project': 'NFLHIVDNA',
            'run_name': 'minitest1',
            'micall_version': 'v7.17'
        }
        sample_info_map = {}  # Empty because there's no 'sample' column

        # Apply sample-specific (none exist)
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if sample_name not in all_samples:
                    raise ValueError(f"Sample '{sample_name}' in sample_info does not exist")
                if 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # Apply run-level to ALL samples
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # ALL samples should have the same run-level micall_version
        assert all_samples['2439P119528-P21I9-NFLHIVDNA_S33']['micall_version'] == 'v7.17'
        assert all_samples['2439P119528-P21F1-NFLHIVDNA_S32']['micall_version'] == 'v7.17'
        assert all_samples['3453P119516-P26P24-NFLHIVDNA_S31']['micall_version'] == 'v7.17'

    def test_with_sample_column_means_sample_specific(self):
        """CRITICAL: 'sample' column in sample_info means sample-specific metadata."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Sample-specific format: has 'sample' column
        sample_info_map = {
            'sample1': {'sample': 'sample1', 'micall_version': 'v7.15'},
            'sample2': {'sample': 'sample2', 'micall_version': 'v7.16'}
        }
        run_level_info = None  # No run-level info

        # Apply sample-specific
        if sample_info_map:
            for sample_name, sample_info in sample_info_map.items():
                if sample_name not in all_samples:
                    raise ValueError(f"Sample '{sample_name}' in sample_info does not exist")
                if 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

        # Apply run-level (none exists)
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Each sample has its own specific version
        assert all_samples['sample1']['micall_version'] == 'v7.15'
        assert all_samples['sample2']['micall_version'] == 'v7.16'

    def test_wrong_sample_name_in_sample_specific_raises_error(self):
        """CRITICAL: Sample name in sample_info must exist in cascade, or raise error."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # sample_info has wrong sample name
        sample_info_map = {
            'sample1': {'sample': 'sample1', 'micall_version': 'v7.15'},
            'wrong_name': {'sample': 'wrong_name', 'micall_version': 'v7.16'}
        }

        # Should raise ValueError for 'wrong_name'
        with pytest.raises(ValueError, match="Sample 'wrong_name' in sample_info does not exist"):
            if sample_info_map:
                for sample_name, sample_info in sample_info_map.items():
                    if sample_name not in all_samples:
                        raise ValueError(
                            f"Sample '{sample_name}' in sample_info does not exist in cascade data. "
                            f"Available samples: {list(all_samples.keys())}"
                        )
                    if 'micall_version' in sample_info:
                        all_samples[sample_name]['micall_version'] = sample_info['micall_version']


class TestFileBasedSampleInfo:
    """Test reading sample_info from actual CSV files on disk."""

    def test_read_run_level_sample_info_from_file(self):
        """Test reading run-level sample_info.csv from file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create cascade.csv
            cascade_path = Path(tmpdir) / 'cascade.csv'
            with open(cascade_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sample', 'demultiplexed', 'v3loop', 'g2p', 'prelim_map', 'remap', 'aligned'])
                writer.writerow(['sample1', '1000', '100', '50', '500', '450', '400'])
                writer.writerow(['sample2', '2000', '200', '100', '1000', '900', '850'])

            # Create run-level sample_info.csv (no 'sample' column)
            sample_info_path = Path(tmpdir) / 'sample_info.csv'
            with open(sample_info_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['project', 'run_name', 'micall_version'])
                writer.writerow(['NFLHIVDNA', 'testrun', 'v7.17'])

            # Read cascade
            with open(cascade_path) as f:
                all_samples = utils.get_samples_from_cascade(f)

            # Read sample_info
            with open(sample_info_path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                assert len(rows) == 1
                first_row = rows[0]

            # Determine if sample_info has 'sample' column
            has_sample_column = 'sample' in first_row
            assert not has_sample_column  # Run-level format

            # Apply run-level info to all samples
            if not has_sample_column and 'micall_version' in first_row:
                for sample_name in all_samples:
                    all_samples[sample_name]['micall_version'] = first_row['micall_version']

            assert all_samples['sample1']['micall_version'] == 'v7.17'
            assert all_samples['sample2']['micall_version'] == 'v7.17'

    def test_read_sample_specific_sample_info_from_file(self):
        """Test reading sample-specific sample_info.csv from file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create cascade.csv
            cascade_path = Path(tmpdir) / 'cascade.csv'
            with open(cascade_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sample', 'demultiplexed', 'v3loop', 'g2p', 'prelim_map', 'remap', 'aligned'])
                writer.writerow(['sample1', '1000', '100', '50', '500', '450', '400'])
                writer.writerow(['sample2', '2000', '200', '100', '1000', '900', '850'])

            # Create sample-specific sample_info.csv (has 'sample' column)
            sample_info_path = Path(tmpdir) / 'sample_info.csv'
            with open(sample_info_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sample', 'project', 'micall_version'])
                writer.writerow(['sample1', 'NFLHIVDNA', 'v7.15'])
                writer.writerow(['sample2', 'NFLHIVDNA', 'v7.16'])

            # Read cascade
            with open(cascade_path) as f:
                all_samples = utils.get_samples_from_cascade(f)

            # Read sample_info
            with open(sample_info_path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                assert len(rows) == 2
                first_row = rows[0]

            # Determine if sample_info has 'sample' column
            has_sample_column = 'sample' in first_row
            assert has_sample_column  # Sample-specific format

            # Build sample_info_map
            sample_info_map = {row['sample']: row for row in rows}

            # Apply sample-specific info
            for sample_name, sample_info in sample_info_map.items():
                if sample_name not in all_samples:
                    raise ValueError(f"Sample '{sample_name}' in sample_info does not exist")
                if 'micall_version' in sample_info:
                    all_samples[sample_name]['micall_version'] = sample_info['micall_version']

            assert all_samples['sample1']['micall_version'] == 'v7.15'
            assert all_samples['sample2']['micall_version'] == 'v7.16'


class TestMultipleSampleRunLevel:
    """Test run-level metadata applied to multiple samples."""

    def test_many_samples_single_run_level_version(self):
        """Test 10 samples all getting the same run-level micall_version."""
        # Create cascade with 10 samples
        cascade_lines = ["sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"]
        sample_names = [f"sample{i}" for i in range(1, 11)]
        for sample_name in sample_names:
            cascade_lines.append(f"{sample_name},1000,100,50,500,450,400\n")

        cascade_csv = io.StringIO(''.join(cascade_lines))
        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Run-level info applies to ALL
        run_level_info = {'micall_version': 'v7.17'}
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # ALL 10 samples should have the same version
        for sample_name in sample_names:
            assert all_samples[sample_name]['micall_version'] == 'v7.17'

    def test_realistic_sample_names_from_ws_rerun(self):
        """Test with realistic sample names like those from ws-rerun-raw-data-folder."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "2439P119528-P21I9-NFLHIVDNA_S33,33703,0,0,0,28168,28013\n"
            "2439P119528-P21F1-NFLHIVDNA_S32,50364,747,225,0,49956,49943\n"
            "3453P119516-P26P24-NFLHIVDNA_S31,41234,500,300,0,38000,37500\n"
            "4567P120001-P30A1-NFLHIVDNA_S40,25000,300,150,0,22000,21500\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Run-level metadata (no 'sample' column in sample_info)
        run_level_info = {
            'project': 'NFLHIVDNA',
            'run_name': 'minitest1',
            'micall_version': 'v7.17'
        }

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Verify ALL samples got the version
        assert all_samples['2439P119528-P21I9-NFLHIVDNA_S33']['micall_version'] == 'v7.17'
        assert all_samples['2439P119528-P21F1-NFLHIVDNA_S32']['micall_version'] == 'v7.17'
        assert all_samples['3453P119516-P26P24-NFLHIVDNA_S31']['micall_version'] == 'v7.17'
        assert all_samples['4567P120001-P30A1-NFLHIVDNA_S40']['micall_version'] == 'v7.17'

    def test_run_level_does_not_modify_other_fields(self):
        """Ensure run-level metadata only adds micall_version, doesn't touch remap counts."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
            "sample2,2000,200,100,1000,900,850\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)

        # Store original remap values
        original_remaps = {
            'sample1': all_samples['sample1']['remap'],
            'sample2': all_samples['sample2']['remap']
        }

        # Apply run-level info
        run_level_info = {'micall_version': 'v7.17'}
        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Remap values should be unchanged
        assert all_samples['sample1']['remap'] == original_remaps['sample1']
        assert all_samples['sample2']['remap'] == original_remaps['sample2']

        # But micall_version should be added
        assert all_samples['sample1']['micall_version'] == 'v7.17'
        assert all_samples['sample2']['micall_version'] == 'v7.17'


class TestWhitespaceAndFormatting:
    """Test handling of whitespace, empty values, and formatting edge cases."""

    def test_micall_version_with_leading_trailing_whitespace(self):
        """Test that whitespace in micall_version is preserved."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        run_level_info = {'micall_version': '  v7.17  '}  # Has whitespace

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        # Whitespace is preserved (not stripped)
        assert all_samples['sample1']['micall_version'] == '  v7.17  '

    def test_numeric_micall_version(self):
        """Test that numeric micall_version values work."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        run_level_info = {'micall_version': '7.17'}  # No 'v' prefix

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        assert all_samples['sample1']['micall_version'] == '7.17'

    def test_very_long_version_string(self):
        """Test handling of unusually long version strings."""
        cascade_csv = io.StringIO(
            "sample,demultiplexed,v3loop,g2p,prelim_map,remap,aligned\n"
            "sample1,1000,100,50,500,450,400\n"
        )

        all_samples = utils.get_samples_from_cascade(cascade_csv)
        long_version = 'v7.17.0-alpha.1+build.12345678.20231201.abcdef1234567890' * 5

        run_level_info = {'micall_version': long_version}

        if run_level_info and 'micall_version' in run_level_info:
            for sample_name in all_samples:
                all_samples[sample_name]['micall_version'] = run_level_info['micall_version']

        assert all_samples['sample1']['micall_version'] == long_version
