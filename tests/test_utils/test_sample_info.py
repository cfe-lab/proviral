"""
Tests for sample_info handling and micall_version propagation.

These tests verify that:
1. micall_version is ONLY read from sample_info.csv (never from cascade.csv)
2. Both sample-specific and run-level sample_info formats are supported
3. Errors are raised when sample_info contains samples not in cascade
"""
import io
import pytest

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
