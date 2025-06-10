"""
Test suite for the compare_runs.py script.

This module contains comprehensive tests for the proviral pipeline comparison
functionality, including:
- Basic file comparison logic
- Directory structure validation
- Missing file handling
- CSV content comparison
- Location field standardization
- JSON output formatting
- Error handling scenarios
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch

from scripts.compare_runs.__main__ import (
    find_versions,
    get_proviral_csv_files,
    read_csv_file,
    compare_csv_contents,
    compare_proviral_files,
    compare_runs,
    main,
    _determine_row_difference_severity,
    _determine_row_difference_confidence,
)
from scripts.compare_runs.discrepancy import (
    Severity,
    Confidence,
    ComparisonReport,
    RowDifferenceDiscrepancy,
    MissingFileDiscrepancy,
    MissingDirectoryDiscrepancy,
    HeaderDifferenceDiscrepancy,
    RowCountDifferenceDiscrepancy,
    ColumnCountDifferenceDiscrepancy,
    FileReadErrorDiscrepancy,
    NoIndexColumnDiscrepancy,
    DuplicateColumnNamesDiscrepancy,
    ColumnOrderDifferenceDiscrepancy,
    MissingRowDiscrepancy,
    ExtraRowDiscrepancy,
)
from scripts.compare_runs.index_discovery import (
    find_unique_value_columns,
    count_shared_values,
    _count_shared_values,
    discover_index_column,
)


class TestDiscrepancy:
    """Test the Discrepancy class functionality."""

    def test_discrepancy_creation(self):
        """Test basic discrepancy creation."""
        # Changed to use a specific discrepancy type: RowDifferenceDiscrepancy
        discrepancy = RowDifferenceDiscrepancy(
            severity=Severity.HIGH,
            confidence=Confidence.HIGH,
            description="Test discrepancy",
            file="test.csv",
            version="version_7.15",
            run="both",  # From old location["run"]
            row=2,  # From old location["row"]
            # Provide placeholder values for new mandatory fields of RowDifferenceDiscrepancy
            index_column="placeholder_idx_col",
            index_value="placeholder_idx_val",
            position_run1=1,
            position_run2=1,
            changed_columns=["placeholder_col"],
            total_field_changes=1,
            field_changes={
                "placeholder_col": {"run1": "A", "run2": "B", "change_type": "text"}
            },
            change_types=["text"],
        )

        assert isinstance(discrepancy, RowDifferenceDiscrepancy)
        assert discrepancy.severity == Severity.HIGH
        assert discrepancy.confidence == Confidence.HIGH
        assert discrepancy.description == "Test discrepancy"
        # Assertions for location and values properties (should still work due to base class implementation)
        assert discrepancy.location["file"] == "test.csv"
        assert discrepancy.location["row"] == 2

    def test_discrepancy_to_dict(self):
        """Test conversion of discrepancy to dictionary."""
        # Changed to use a specific discrepancy type: MissingFileDiscrepancy
        discrepancy = MissingFileDiscrepancy(
            severity=Severity.CRITICAL,
            confidence=Confidence.HIGH,
            description="File missing",
            file="missing.csv",  # From old location
            version="version_7.15",  # From old location
            run="run1",  # From old location
            missing_from="run1",  # From old values
            present_in="run2",  # Inferred: if missing from run1, present in run2
        )

        result = discrepancy.to_dict()

        # Type field removed; ensure severity and confidence only
        assert result["severity"] == "CRITICAL"
        assert result["confidence"] == "HIGH"
        assert result["description"] == "File missing"
        assert result["location"]["file"] == "missing.csv"
        assert (
            result["location"]["missing_from"] == "run1"
        )  # Added by MissingFileDiscrepancy._add_location_fields
        assert (
            result["missing_from"] == "run1"
        )  # Added by MissingFileDiscrepancy._add_values_fields
        assert (
            result["present_in"] == "run2"
        )  # Added by MissingFileDiscrepancy._add_values_fields


class TestComparisonReport:
    """Test the ComparisonReport class functionality."""

    def test_report_initialization(self):
        """Test basic report initialization."""
        run1_dir = Path("/fake/run1")
        run2_dir = Path("/fake/run2")
        report = ComparisonReport(run1_dir, run2_dir)

        assert report.run1_dir == run1_dir
        assert report.run2_dir == run2_dir
        assert report.versions_run1 == []
        assert report.versions_run2 == []
        assert report.common_versions == []
        assert isinstance(report.timestamp, str)

    def test_add_discrepancy(self):
        """Test adding discrepancies to report."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))
        # Changed to use a specific discrepancy type: RowDifferenceDiscrepancy
        discrepancy = RowDifferenceDiscrepancy(
            severity=Severity.MEDIUM,
            confidence=Confidence.HIGH,
            description="Row differs",
            file="test.csv",
            version="version_7.15",
            run="both",  # From old location["run"]
            row=2,  # From old location["row"]
            # Provide placeholder values for new mandatory fields
            index_column="placeholder_idx_col",
            index_value="placeholder_idx_val",
            position_run1=1,
            position_run2=1,
            changed_columns=[],
            total_field_changes=0,
            field_changes={},
            change_types=[],
            # No _original_values needed as the old test didn't provide a 'values' dict
        )

        report.add_discrepancy("version_7.15", "test.csv", discrepancy)

        assert "version_7.15" in report.results
        assert "test.csv" in report.results["version_7.15"]
        assert report.results["version_7.15"]["test.csv"]["status"] == "differs"
        assert len(report.results["version_7.15"]["test.csv"]["discrepancies"]) == 1

    def test_mark_file_identical(self):
        """Test marking files as identical."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))

        report.mark_file_identical("version_7.15", "identical.csv")

        assert report.results["version_7.15"]["identical.csv"]["status"] == "identical"
        assert report.results["version_7.15"]["identical.csv"]["discrepancies"] == []

    def test_get_summary(self):
        """Test summary statistics generation."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))

        # Add various discrepancies using specific types
        critical_discrepancy = MissingFileDiscrepancy(
            severity=Severity.CRITICAL,
            confidence=Confidence.HIGH,
            description="Critical",
            file="file1.csv",
            version="version_7.15",
            run="run1",  # run where it's missing
            missing_from="run1",
            present_in="run2",  # Placeholders
        )
        high_discrepancy = HeaderDifferenceDiscrepancy(
            severity=Severity.HIGH,
            confidence=Confidence.HIGH,
            description="High",
            file="file2.csv",
            version="version_7.15",
            run="both",  # Comparison of both
            # Placeholders for HeaderDifferenceDiscrepancy fields
            changed_headers=[],
            total_header_changes=0,
            header_changes={},
            run1_header_count=0,
            run2_header_count=0,
        )
        medium_discrepancy = RowDifferenceDiscrepancy(
            severity=Severity.MEDIUM,
            confidence=Confidence.MEDIUM,
            description="Medium",
            file="file3.csv",
            version="version_7.15",
            run="both",  # Comparison of both
            row=1,  # Placeholder
            # Placeholders for RowDifferenceDiscrepancy fields
            index_column="idx",
            index_value="val",
            position_run1=1,
            position_run2=1,
            changed_columns=[],
            total_field_changes=0,
            field_changes={},
            change_types=[],
        )

        report.add_discrepancy("version_7.15", "file1.csv", critical_discrepancy)
        report.add_discrepancy("version_7.15", "file2.csv", high_discrepancy)
        report.add_discrepancy("version_7.15", "file3.csv", medium_discrepancy)

        summary = report.get_summary()

        assert summary["total_discrepancies"] == 3
        assert summary["critical_discrepancies"] == 1
        assert summary["high_discrepancies"] == 1
        assert summary["medium_discrepancies"] == 1
        assert summary["low_discrepancies"] == 0

    def test_to_dict(self):
        """Test conversion to dictionary for JSON output."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))
        report.versions_run1 = ["version_7.15"]
        report.versions_run2 = ["version_7.15", "version_7.17"]
        report.common_versions = ["version_7.15"]

        result = report.to_dict()

        assert "metadata" in result
        assert "summary" in result
        assert "results" in result
        assert result["metadata"]["run1_dir"] == "/fake/run1"
        assert result["metadata"]["versions_run1"] == ["version_7.15"]

    def test_to_json(self):
        """Test JSON serialization."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))

        json_str = report.to_json()
        parsed = json.loads(json_str)

        assert "metadata" in parsed
        assert "summary" in parsed
        assert "results" in parsed

    def test_set_file_index_info(self):
        """Test setting index column information for a file."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))
        index_info = {
            "index_column": 0,
            "column_name": "sample_id",
            "shared_values": 5,
            "reason": "success",
        }

        report.set_file_index_info("version_7.15", "test.csv", index_info)

        assert "version_7.15" in report.results
        assert "test.csv" in report.results["version_7.15"]
        assert report.results["version_7.15"]["test.csv"]["index_column"] == index_info

    def test_mark_file_identical_with_index_info(self):
        """Test marking files as identical with index column information."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))
        index_info = {
            "index_column": 1,
            "column_name": "unique_id",
            "shared_values": 10,
            "reason": "success",
        }

        report.mark_file_identical("version_7.15", "identical.csv", index_info)

        assert report.results["version_7.15"]["identical.csv"]["status"] == "identical"
        assert report.results["version_7.15"]["identical.csv"]["discrepancies"] == []
        assert (
            report.results["version_7.15"]["identical.csv"]["index_column"]
            == index_info
        )


class TestFileOperations:
    """Test file and directory operation functions."""

    def test_find_versions_missing_results_dir(self, tmp_path):
        """Test finding versions when Results directory doesn't exist."""
        run_dir = tmp_path / "fake_run"
        run_dir.mkdir()

        versions = find_versions(run_dir)

        assert versions == []

    def test_find_versions_success(self, tmp_path):
        """Test successful version discovery."""
        run_dir = tmp_path / "fake_run"
        results_dir = run_dir / "Results"
        results_dir.mkdir(parents=True)

        # Create version directories
        (results_dir / "version_7.15").mkdir()
        (results_dir / "version_7.17").mkdir()
        (results_dir / "not_a_version").mkdir()  # Should be ignored

        versions = find_versions(run_dir)

        assert versions == ["version_7.15", "version_7.17"]

    def test_get_proviral_csv_files_missing_dir(self, tmp_path):
        """Test getting CSV files when proviral directory doesn't exist."""
        results_dir = tmp_path / "Results"
        results_dir.mkdir()

        csv_files = get_proviral_csv_files(results_dir, "version_7.15")

        assert csv_files == {}

    def test_get_proviral_csv_files_success(self, tmp_path):
        """Test successful CSV file discovery."""
        results_dir = tmp_path / "Results"
        proviral_dir = results_dir / "version_7.15" / "proviral"
        proviral_dir.mkdir(parents=True)

        # Create CSV files
        (proviral_dir / "test1.csv").touch()
        (proviral_dir / "test2.csv").touch()
        (proviral_dir / "not_csv.txt").touch()  # Should be ignored

        csv_files = get_proviral_csv_files(results_dir, "version_7.15")

        assert len(csv_files) == 2
        assert "test1.csv" in csv_files
        assert "test2.csv" in csv_files
        assert "not_csv.txt" not in csv_files

    def test_read_csv_file_success(self, tmp_path):
        """Test successful CSV file reading."""
        csv_file = tmp_path / "test.csv"
        csv_file.write_text("header1,header2\nvalue1,value2\nvalue3,value4\n")

        content = read_csv_file(csv_file)

        assert len(content) == 3
        assert content[0] == ["header1", "header2"]
        assert content[1] == ["value1", "value2"]
        assert content[2] == ["value3", "value4"]

    def test_read_csv_file_nonexistent(self, tmp_path):
        """Test reading non-existent CSV file."""
        csv_file = tmp_path / "nonexistent.csv"

        content = read_csv_file(csv_file)

        assert content == []


class TestCSVComparison:
    """Test CSV content comparison functions."""

    def test_compare_identical_files(self, tmp_path):
        """Test comparison of identical CSV files."""
        # Create identical files
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        content = "header1,header2\nvalue1,value2\n"
        file1.write_text(content)
        file2.write_text(content)

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        assert discrepancies == []

    def test_compare_different_row_count(self, tmp_path):
        """Test comparison with different row counts."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header2\nvalue1,value2\nvalue3,value4\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find row count difference
        row_count_discrepancy = next(
            (d for d in discrepancies if isinstance(d, RowCountDifferenceDiscrepancy)),
            None,
        )
        assert row_count_discrepancy is not None
        assert row_count_discrepancy.severity == Severity.HIGH
        assert row_count_discrepancy.location["file"] == "test.csv"
        assert row_count_discrepancy.location["version"] == "version_7.15"
        assert row_count_discrepancy.location["run"] == "both"

    def test_compare_different_headers(self, tmp_path):
        """Test comparison with different headers."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header3\nvalue1,value2\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find header difference
        header_discrepancy = next(
            (d for d in discrepancies if isinstance(d, HeaderDifferenceDiscrepancy)),
            None,
        )
        assert header_discrepancy is not None
        assert header_discrepancy.severity == Severity.CRITICAL
        assert header_discrepancy.location["file"] == "test.csv"
        assert header_discrepancy.location["row"] == 1

    def test_compare_different_data_rows(self, tmp_path):
        """Test comparison with different data rows."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header2\nvalue1,different\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find row difference
        row_discrepancy = next(
            (d for d in discrepancies if isinstance(d, RowDifferenceDiscrepancy)), None
        )
        assert row_discrepancy is not None
        assert row_discrepancy.location["file"] == "test.csv"
        assert row_discrepancy.location["row"] == 2

    def test_compare_empty_files(self, tmp_path):
        """Test comparison of empty files."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("")
        file2.write_text("")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        assert discrepancies == []

    def test_compare_one_empty_file(self, tmp_path):
        """Test comparison where one file is empty."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\n")
        file2.write_text("")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find file read error
        error_discrepancy = next(
            (d for d in discrepancies if isinstance(d, FileReadErrorDiscrepancy)),
            None,
        )
        assert error_discrepancy is not None
        assert error_discrepancy.severity == Severity.CRITICAL
        assert error_discrepancy.location["run"] == "run2"


class TestSeverityDetermination:
    """Test severity and confidence determination functions."""

    def test_determine_row_difference_severity_critical_outcome(self):
        """Test severity determination for critical outcome differences."""
        row1 = ["sample1", "success", "0.95"]
        row2 = ["sample1", "failed", "0.95"]

        # Mock column differences that indicate an outcome change
        column_differences = {
            "indices": [1],
            "change_types": ["outcome"],
            "field_changes": [{"change_type": "outcome"}],
        }

        severity = _determine_row_difference_severity(row1, row2, 1, column_differences)

        assert severity == Severity.HIGH

    def test_determine_row_difference_severity_numeric_difference(self):
        """Test severity determination for numeric differences."""
        row1 = ["sample1", "success", "0.95"]
        row2 = ["sample1", "success", "0.85"]

        # Mock column differences that indicate a numeric change
        column_differences = {
            "indices": [2],
            "change_types": ["numeric"],
            "field_changes": [{"change_type": "numeric"}],
        }

        severity = _determine_row_difference_severity(row1, row2, 1, column_differences)

        assert severity == Severity.MEDIUM

    def test_determine_row_difference_severity_minor_difference(self):
        """Test severity determination for minor differences."""
        row1 = ["sample1", "other", "text"]
        row2 = ["sample1", "different", "text"]

        # Mock column differences that indicate a text change
        column_differences = {
            "indices": [1],
            "change_types": ["text"],
            "field_changes": [{"change_type": "text"}],
        }

        severity = _determine_row_difference_severity(row1, row2, 1, column_differences)

        assert severity == Severity.MEDIUM

    def test_determine_row_difference_confidence_high(self):
        """Test confidence determination for clear differences."""
        row1 = ["a", "b", "c"]
        row2 = ["x", "y", "z"]

        confidence = _determine_row_difference_confidence(row1, row2)

        assert confidence == Confidence.HIGH

    def test_determine_row_difference_confidence_medium(self):
        """Test confidence determination for single differences."""
        row1 = ["a", "b", "c"]
        row2 = ["a", "b", "z"]

        confidence = _determine_row_difference_confidence(row1, row2)

        assert confidence == Confidence.MEDIUM


class TestProviralFileComparison:
    """Test the compare_proviral_files function."""

    def test_compare_missing_directories(self, tmp_path):
        """Test comparison when proviral directories are missing."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        run1_dir.mkdir()
        run2_dir.mkdir()

        report = ComparisonReport(run1_dir, run2_dir)

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        # Should have missing directory discrepancy
        assert "version_7.15" in report.results
        assert "proviral_directory" in report.results["version_7.15"]

    def test_compare_missing_file_in_run1(self, tmp_path):
        """Test comparison when a file is missing in run1."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        proviral1 = run1_dir / "Results" / "version_7.15" / "proviral"
        proviral2 = run2_dir / "Results" / "version_7.15" / "proviral"
        proviral1.mkdir(parents=True)
        proviral2.mkdir(parents=True)

        # Create file only in run2
        (proviral2 / "missing.csv").write_text("header\nvalue\n")

        report = ComparisonReport(run1_dir, run2_dir)

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        # Should have missing file discrepancy
        assert "missing.csv" in report.results["version_7.15"]
        discrepancy = report.results["version_7.15"]["missing.csv"]["discrepancies"][0]
        assert discrepancy["location"]["missing_from"] == "run1"

    def test_compare_identical_files(self, tmp_path):
        """Test comparison of identical files."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        proviral1 = run1_dir / "Results" / "version_7.15" / "proviral"
        proviral2 = run2_dir / "Results" / "version_7.15" / "proviral"
        proviral1.mkdir(parents=True)
        proviral2.mkdir(parents=True)

        # Create identical files
        content = "header1,header2\nvalue1,value2\n"
        (proviral1 / "identical.csv").write_text(content)
        (proviral2 / "identical.csv").write_text(content)

        report = ComparisonReport(run1_dir, run2_dir)

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        # Should mark file as identical
        assert report.results["version_7.15"]["identical.csv"]["status"] == "identical"


class TestFullComparison:
    """Test the complete compare_runs function."""

    def test_compare_no_common_versions(self, tmp_path):
        """Test comparison when there are no common versions."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        results1 = run1_dir / "Results"
        results2 = run2_dir / "Results"
        results1.mkdir(parents=True)
        results2.mkdir(parents=True)

        # Create different versions
        (results1 / "version_7.15").mkdir()
        (results2 / "version_7.17").mkdir()

        report = compare_runs(run1_dir, run2_dir)

        # Should have no common versions discrepancy
        assert "metadata" in report.results
        assert "versions" in report.results["metadata"]

    def test_compare_version_mismatches(self, tmp_path):
        """Test comparison with version mismatches."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        results1 = run1_dir / "Results"
        results2 = run2_dir / "Results"
        results1.mkdir(parents=True)
        results2.mkdir(parents=True)

        # Create versions with some overlap
        (results1 / "version_7.15").mkdir()
        (results1 / "version_7.16").mkdir()
        (results2 / "version_7.15").mkdir()
        (results2 / "version_7.17").mkdir()

        report = compare_runs(run1_dir, run2_dir)

        assert report.common_versions == ["version_7.15"]
        assert "metadata" in report.results

    def test_compare_successful_run(self, tmp_path):
        """Test a successful comparison run."""
        # Setup full directory structure
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"

        for run_dir in [run1_dir, run2_dir]:
            proviral_dir = run_dir / "Results" / "version_7.15" / "proviral"
            proviral_dir.mkdir(parents=True)
            (proviral_dir / "test.csv").write_text("header\nvalue\n")

        # Make files different
        (run2_dir / "Results" / "version_7.15" / "proviral" / "test.csv").write_text(
            "header\ndifferent\n"
        )

        report = compare_runs(run1_dir, run2_dir)

        assert report.common_versions == ["version_7.15"]
        assert "version_7.15" in report.results
        assert "test.csv" in report.results["version_7.15"]

    def test_compare_with_comparison_error(self, tmp_path):
        """Test handling of comparison errors."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        results1 = run1_dir / "Results"
        results2 = run2_dir / "Results"
        results1.mkdir(parents=True)
        results2.mkdir(parents=True)
        (results1 / "version_7.15").mkdir()
        (results2 / "version_7.15").mkdir()

        # Create proviral directories with files to ensure compare_proviral_files gets called
        proviral1 = results1 / "version_7.15" / "proviral"
        proviral2 = results2 / "version_7.15" / "proviral"
        proviral1.mkdir()
        proviral2.mkdir()
        (proviral1 / "test.csv").write_text("header\nvalue\n")
        (proviral2 / "test.csv").write_text("header\nvalue\n")

        # Mock compare_proviral_files to raise an exception - use the correct module path
        with patch(
            "scripts.compare_runs.__main__.compare_proviral_files",
            side_effect=Exception("Test error"),
        ):
            report = compare_runs(run1_dir, run2_dir)

        # Should have comparison error discrepancy
        assert "comparison_error" in report.results["version_7.15"]
        discrepancy = report.results["version_7.15"]["comparison_error"][
            "discrepancies"
        ][0]

        assert "Test error" in discrepancy["description"]


class TestLocationFieldStandardization:
    """Test that location fields are standardized across all discrepancy types."""

    def test_missing_file_location_fields(self, tmp_path):
        """Test location fields for missing file discrepancies."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        proviral1 = run1_dir / "Results" / "version_7.15" / "proviral"
        proviral2 = run2_dir / "Results" / "version_7.15" / "proviral"
        proviral1.mkdir(parents=True)
        proviral2.mkdir(parents=True)

        # Create file only in run2
        (proviral2 / "missing.csv").write_text("header\nvalue\n")

        report = ComparisonReport(run1_dir, run2_dir)
        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        discrepancy = report.results["version_7.15"]["missing.csv"]["discrepancies"][0]
        location = discrepancy["location"]

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "missing_from" in location  # Type-specific field

        assert location["file"] == "missing.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "run1"

    def test_row_difference_location_fields(self, tmp_path):
        """Test that NO_INDEX_COLUMN discrepancies are reported when no suitable index column is available."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue1\n")
        file2.write_text("header\nvalue2\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find NO_INDEX_COLUMN discrepancy instead of row differences
        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )
        assert no_index_discrepancy is not None

        location = no_index_discrepancy.location

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "reason" in location  # Type-specific field

        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"

    def test_file_read_error_location_fields(self, tmp_path):
        """Test location fields for file read error discrepancies."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue\n")
        file2.write_text("")  # Empty file

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")
        error_discrepancy = next(
            (d for d in discrepancies if isinstance(d, FileReadErrorDiscrepancy)),
            None,
        )

        location = error_discrepancy.location

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "file_path" in location  # Type-specific field

        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "run2"


class TestMainFunction:
    """Test the main function and command-line interface."""

    def test_main_with_nonexistent_directory(self, tmp_path, capsys):
        """Test main function with non-existent directory."""
        nonexistent = tmp_path / "nonexistent"
        existing = tmp_path / "existing"
        existing.mkdir()

        with patch("sys.argv", ["compare_runs", str(nonexistent), str(existing)]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1

    def test_main_successful_comparison(self, tmp_path, capsys):
        """Test successful main function execution."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"

        for run_dir in [run1_dir, run2_dir]:
            proviral_dir = run_dir / "Results" / "version_7.15" / "proviral"
            proviral_dir.mkdir(parents=True)
            (proviral_dir / "test.csv").write_text("header\nvalue\n")

        with patch("sys.argv", ["compare_runs", str(run1_dir), str(run2_dir)]):
            main()

        captured = capsys.readouterr()
        output = json.loads(captured.out)

        assert "metadata" in output
        assert "summary" in output
        assert "results" in output

    def test_main_with_output_file(self, tmp_path):
        """Test main function with output file option."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        output_file = tmp_path / "output.json"

        for run_dir in [run1_dir, run2_dir]:
            proviral_dir = run_dir / "Results" / "version_7.15" / "proviral"
            proviral_dir.mkdir(parents=True)
            (proviral_dir / "test.csv").write_text("header\nvalue\n")

        with patch(
            "sys.argv",
            ["compare_runs", str(run1_dir), str(run2_dir), "-o", str(output_file)],
        ):
            main()

        assert output_file.exists()
        output = json.loads(output_file.read_text())
        assert "metadata" in output

    def test_main_with_compact_output(self, tmp_path, capsys):
        """Test main function with compact JSON output."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"

        for run_dir in [run1_dir, run2_dir]:
            proviral_dir = run_dir / "Results" / "version_7.15" / "proviral"
            proviral_dir.mkdir(parents=True)
            (proviral_dir / "test.csv").write_text("header\nvalue\n")

        with patch(
            "sys.argv", ["compare_runs", str(run1_dir), str(run2_dir), "--compact"]
        ):
            main()

        captured = capsys.readouterr()
        # Compact JSON should not have newlines except at the end
        assert captured.out.count("\n") == 1

    def test_main_keyboard_interrupt(self, tmp_path):
        """Test main function handling keyboard interrupt."""
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        run1_dir.mkdir()
        run2_dir.mkdir()

        # Create Results directories to ensure compare_runs gets called
        (run1_dir / "Results").mkdir()
        (run2_dir / "Results").mkdir()

        with patch(
            "scripts.compare_runs.__main__.compare_runs",
            side_effect=KeyboardInterrupt(),
        ):
            with patch("sys.argv", ["compare_runs", str(run1_dir), str(run2_dir)]):
                with pytest.raises(SystemExit) as exc_info:
                    main()
                assert exc_info.value.code == 1


class TestIndexColumnDiscovery:
    """Test index column discovery functionality."""

    def test_find_unique_value_columns_empty_data(self):
        """Test find_unique_value_columns with empty data."""
        assert find_unique_value_columns([]) == []
        assert find_unique_value_columns([["header"]]) == []  # Only header, no data

    def test_find_unique_value_columns_all_unique(self):
        """Test find_unique_value_columns with all unique columns."""
        csv_data = [
            ["id", "name", "value"],
            ["1", "alice", "100"],
            ["2", "bob", "200"],
            ["3", "charlie", "300"],
        ]
        unique_columns = find_unique_value_columns(csv_data)
        # All columns should be unique
        assert sorted(unique_columns) == ["id", "name", "value"]

    def test_find_unique_value_columns_with_duplicates(self):
        """Test find_unique_value_columns with some duplicate values."""
        csv_data = [
            ["id", "category", "value"],
            ["1", "A", "100"],
            ["2", "A", "200"],  # category "A" is duplicate
            ["3", "B", "300"],
        ]
        unique_columns = find_unique_value_columns(csv_data)
        # Only id and value columns should be unique
        assert sorted(unique_columns) == ["id", "value"]

    def test_find_unique_value_columns_with_empty_values(self):
        """Test find_unique_value_columns ignoring empty values."""
        csv_data = [
            ["id", "optional", "value"],
            ["1", "", "100"],
            ["2", "data", "200"],
            ["3", "", "300"],
        ]
        unique_columns = find_unique_value_columns(csv_data)
        # All columns should be unique (empty values are ignored)
        assert sorted(unique_columns) == ["id", "optional", "value"]

    def test_count_shared_values_no_shared(self):
        """Test count_shared_values with no shared values."""
        csv_data1 = [["id", "value"], ["1", "A"], ["2", "B"]]
        csv_data2 = [["id", "value"], ["3", "C"], ["4", "D"]]
        shared_count = count_shared_values(csv_data1, csv_data2, "id")
        assert shared_count == 0

    def test_count_shared_values_some_shared(self):
        """Test count_shared_values with some shared values."""
        csv_data1 = [["id", "value"], ["1", "A"], ["2", "B"], ["3", "C"]]
        csv_data2 = [
            ["id", "value"],
            ["2", "X"],  # id "2" is shared
            ["3", "Y"],  # id "3" is shared
            ["4", "Z"],
        ]
        shared_count = count_shared_values(csv_data1, csv_data2, "id")
        assert shared_count == 2

    def test_count_shared_values_with_empty_values(self):
        """Test count_shared_values ignoring empty values."""
        csv_data1 = [
            ["id", "value"],
            ["1", "A"],
            ["", "B"],  # Empty id should be ignored
            ["3", "C"],
        ]
        csv_data2 = [
            ["id", "value"],
            ["1", "X"],  # id "1" is shared
            ["", "Y"],  # Empty id should be ignored
            ["4", "Z"],
        ]
        shared_count = _count_shared_values(csv_data1, csv_data2, "id")
        assert shared_count == 1

    def test_discover_index_column_empty_data(self):
        """Test discover_index_column with empty data."""
        assert discover_index_column([], []) is None
        assert discover_index_column([["header"]], []) is None
        assert discover_index_column([], [["header"]]) is None

    def test_discover_index_column_no_unique_columns(self):
        """Test discover_index_column when no columns are unique in both files."""
        csv_data1 = [
            ["id", "category"],
            ["1", "A"],
            ["1", "A"],  # id is duplicate
        ]
        csv_data2 = [
            ["id", "category"],
            ["2", "B"],
            ["2", "B"],  # id is duplicate
        ]
        result = discover_index_column(csv_data1, csv_data2)
        assert result is None

    def test_discover_index_column_no_shared_values(self):
        """Test discover_index_column when unique columns exist but no shared values."""
        csv_data1 = [["id", "value"], ["1", "A"], ["2", "B"]]
        csv_data2 = [
            ["id", "value"],
            ["3", "C"],  # Different ids, no shared values
            ["4", "D"],
        ]
        result = discover_index_column(csv_data1, csv_data2)

        assert result is not None
        assert result["index_column"] is None
        assert result["shared_values"] == 0
        assert result["reason"] == "no_shared_values"
        assert result["candidate_columns"] == 2  # Both id and value are unique

    def test_discover_index_column_successful_discovery(self):
        """Test discover_index_column with successful index column discovery."""
        csv_data1 = [
            ["sample_id", "batch", "result"],
            ["S001", "B1", "pass"],
            ["S002", "B1", "fail"],
            ["S003", "B2", "pass"],
        ]
        csv_data2 = [
            ["sample_id", "batch", "result"],
            ["S001", "B1", "pass"],  # S001 shared
            ["S002", "B1", "pass"],  # S002 shared (different result)
            ["S004", "B2", "fail"],  # S004 not shared
        ]
        result = discover_index_column(csv_data1, csv_data2)

        assert result is not None
        assert result["index_column"] == 0  # sample_id column
        assert result["column_name"] == "sample_id"
        assert result["shared_values"] == 2  # S001 and S002
        assert result["reason"] == "success"
        assert "column_analysis" in result

    def test_discover_index_column_chooses_best_column(self):
        """Test discover_index_column chooses column with most shared values."""
        csv_data1 = [
            ["id1", "id2", "value"],
            ["A1", "X1", "100"],
            ["A2", "X2", "200"],
            ["A3", "X3", "300"],
        ]
        csv_data2 = [
            ["id1", "id2", "value"],
            ["A1", "Y1", "150"],  # id1 shared: A1, id2 not shared
            ["A2", "Y2", "250"],  # id1 shared: A2, id2 not shared
            ["B3", "X3", "350"],  # id1 not shared, id2 shared: X3
        ]
        result = discover_index_column(csv_data1, csv_data2)

        assert result is not None
        assert result["index_column"] == 0  # id1 has more shared values (2 vs 1)
        assert result["column_name"] == "id1"
        assert result["shared_values"] == 2

    def test_discover_index_column_with_mismatched_headers(self):
        """Test discover_index_column with different headers between files."""
        csv_data1 = [["sample_id", "result"], ["S001", "pass"], ["S002", "fail"]]
        csv_data2 = [
            ["id", "outcome"],  # Different header names
            ["S001", "unique_value"],  # Make result column have fewer shared values
            ["S003", "another_value"],
        ]
        result = discover_index_column(csv_data1, csv_data2)

        # When headers don't match, we cannot safely determine column correspondence
        # by name, so no index column should be discovered
        assert result is None


class TestCSVComparisonWithIndexDiscovery:
    """Test CSV comparison functionality with index column discovery."""

    def test_compare_csv_with_index_discovery_identical_files(self, tmp_path):
        """Test CSV comparison with index discovery for identical files."""
        content = "sample_id,result\nS001,pass\nS002,fail\n"
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text(content)
        file2.write_text(content)

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should be no discrepancies for identical files
        assert discrepancies == []

    def test_compare_csv_with_no_index_column_discrepancy(self, tmp_path):
        """Test CSV comparison generates NO_INDEX_COLUMN discrepancy when appropriate."""
        # Create files with unique columns but no shared values
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,value\n1,A\n2,B\n")
        file2.write_text("id,value\n3,C\n4,D\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find NO_INDEX_COLUMN discrepancy
        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy.severity == Severity.CRITICAL
        assert no_index_discrepancy.confidence == Confidence.HIGH
        assert "Row comparison skipped" in no_index_discrepancy.description

    def test_compare_csv_no_index_discrepancy_when_no_unique_columns(self, tmp_path):
        """Test that NO_INDEX_COLUMN discrepancy is not generated when no unique columns exist."""
        # Create files with no unique columns - all columns have duplicate values
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text(
            "category,status\nA,active\nA,active\n"
        )  # Both columns have duplicates
        file2.write_text(
            "category,status\nB,pending\nB,pending\n"
        )  # Both columns have duplicates

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find NO_INDEX_COLUMN discrepancy since no suitable index column is available
        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy.severity == Severity.CRITICAL
        assert "Row comparison skipped" in no_index_discrepancy.description

    def test_compare_csv_with_successful_index_discovery(self, tmp_path):
        """Test CSV comparison with successful index column discovery."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("sample_id,result\nS001,pass\nS002,fail\n")
        file2.write_text(
            "sample_id,result\nS001,pass\nS002,pass\n"  # S002 result differs
        )

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find row difference but no NO_INDEX_COLUMN discrepancy
        row_discrepancy = next(
            (d for d in discrepancies if isinstance(d, RowDifferenceDiscrepancy)),
            None,
        )
        assert row_discrepancy is not None

        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )
        assert (
            no_index_discrepancy is None
        )  # Should not have this since index was found


class TestIndexColumnLocationFields:
    """Test location fields for NO_INDEX_COLUMN discrepancy type."""

    def test_no_index_column_location_fields(self, tmp_path):
        """Test location fields for NO_INDEX_COLUMN discrepancies."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,value\n1,A\n2,B\n")
        file2.write_text("id,value\n3,C\n4,D\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")
        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )

        assert no_index_discrepancy is not None
        location = no_index_discrepancy.location

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location

        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"


class TestIntegrationWithIndexDiscovery:
    """Test integration of index column discovery with file comparison functions."""

    def test_compare_proviral_files_with_index_info(self, tmp_path):
        """Test that compare_proviral_files stores index column information."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        proviral1 = run1_dir / "Results" / "version_7.15" / "proviral"
        proviral2 = run2_dir / "Results" / "version_7.15" / "proviral"
        proviral1.mkdir(parents=True)
        proviral2.mkdir(parents=True)

        # Create CSV files with index-friendly data
        content1 = "sample_id,result\nS001,pass\nS002,fail\n"
        content2 = (
            "sample_id,result\nS001,pass\nS002,pass\n"  # Different result for S002
        )
        (proviral1 / "test.csv").write_text(content1)
        (proviral2 / "test.csv").write_text(content2)

        report = ComparisonReport(run1_dir, run2_dir)
        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        # Check that index column information is stored
        file_result = report.results["version_7.15"]["test.csv"]
        assert "index_column" in file_result
        index_info = file_result["index_column"]

        assert index_info is not None
        assert index_info["index_column"] == 0  # sample_id column
        assert index_info["column_name"] == "sample_id"
        assert index_info["shared_values"] == 2
        assert index_info["reason"] == "success"

    def test_compare_identical_files_with_index_info(self, tmp_path):
        """Test that identical files also get index column information."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"
        proviral1 = run1_dir / "Results" / "version_7.15" / "proviral"
        proviral2 = run2_dir / "Results" / "version_7.15" / "proviral"
        proviral1.mkdir(parents=True)
        proviral2.mkdir(parents=True)

        # Create identical CSV files
        content = "sample_id,result\nS001,pass\nS002,fail\n"
        (proviral1 / "identical.csv").write_text(content)
        (proviral2 / "identical.csv").write_text(content)

        report = ComparisonReport(run1_dir, run2_dir)
        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report)

        # Check that file is marked identical and has index column info
        file_result = report.results["version_7.15"]["identical.csv"]
        assert file_result["status"] == "identical"
        assert "index_column" in file_result

        index_info = file_result["index_column"]
        assert index_info is not None
        # Both columns are unique and have same shared values, algorithm picks "result" (sorted order)
        assert index_info["index_column"] == 1
        assert index_info["column_name"] == "result"


class TestColumnValidation:
    """Test column name validation and column order comparison functionality."""

    def test_duplicate_column_names_run1(self, tmp_path):
        """Test detection of duplicate column names in run1."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text(
            "header1,header2,header1\nvalue1,value2,value3\n"
        )  # Duplicate header1
        file2.write_text("header1,header2,header3\nvalue1,value2,value3\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find duplicate column names discrepancy
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, DuplicateColumnNamesDiscrepancy)
            ),
            None,
        )
        assert dup_discrepancy is not None
        assert dup_discrepancy.severity == Severity.CRITICAL
        assert dup_discrepancy.confidence == Confidence.HIGH
        assert "Duplicate column name 'header1'" in dup_discrepancy.description
        assert "run1" in dup_discrepancy.description

        # Check location fields
        location = dup_discrepancy.location
        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "run1"
        assert location["duplicate_column"] == "header1"
        assert location["positions"] == [0, 2]  # header1 appears at positions 0 and 2

        # Check values
        assert dup_discrepancy.duplicate_column == "header1"
        assert dup_discrepancy.all_headers == ["header1", "header2", "header1"]

    def test_duplicate_column_names_run2(self, tmp_path):
        """Test detection of duplicate column names in run2."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2,header3\nvalue1,value2,value3\n")
        file2.write_text(
            "col1,col2,col1,col3\nvalue1,value2,value3,value4\n"
        )  # Duplicate col1

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find duplicate column names discrepancy
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, DuplicateColumnNamesDiscrepancy)
            ),
            None,
        )
        assert dup_discrepancy is not None
        assert "Duplicate column name 'col1'" in dup_discrepancy.description
        assert "run2" in dup_discrepancy.description
        assert dup_discrepancy.location["run"] == "run2"
        assert dup_discrepancy.location["duplicate_column"] == "col1"
        assert dup_discrepancy.location["positions"] == [0, 2]

    def test_duplicate_column_names_both_runs(self, tmp_path):
        """Test detection of duplicate column names in both runs."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,id\nvalue1,alice,value3\n")  # Duplicate id
        file2.write_text("name,id,name\nvalue1,value2,value3\n")  # Duplicate name

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find duplicate column names discrepancies for both runs
        dup_discrepancies = [
            d for d in discrepancies if isinstance(d, DuplicateColumnNamesDiscrepancy)
        ]
        assert len(dup_discrepancies) == 2

        # Check that we have one for each run
        run1_discrepancy = next(
            (d for d in dup_discrepancies if d.location["run"] == "run1"), None
        )
        run2_discrepancy = next(
            (d for d in dup_discrepancies if d.location["run"] == "run2"), None
        )

        assert run1_discrepancy is not None
        assert run2_discrepancy is not None
        assert run1_discrepancy.location["duplicate_column"] == "id"
        assert run2_discrepancy.location["duplicate_column"] == "name"

    def test_column_order_difference(self, tmp_path):
        """Test detection of column order differences."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,value\n1,alice,100\n2,bob,200\n")
        file2.write_text(
            "name,id,value\n1,alice,100\n2,bob,200\n"
        )  # Same columns, different order

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find column order difference
        order_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, ColumnOrderDifferenceDiscrepancy)
            ),
            None,
        )
        assert order_discrepancy is not None
        assert order_discrepancy.severity == Severity.MEDIUM
        assert order_discrepancy.confidence == Confidence.HIGH
        assert "Column order differs" in order_discrepancy.description
        assert "2 columns reordered" in order_discrepancy.description

        # Check location fields
        location = order_discrepancy.location
        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"
        assert "reordered_columns" in location
        assert set(location["reordered_columns"]) == {
            "id",
            "name",
        }  # Both id and name changed positions

        # Check values
        assert order_discrepancy.reordered_count == 2

        assert order_discrepancy.order_differences["headers_run1"] == [
            "id",
            "name",
            "value",
        ]
        assert order_discrepancy.order_differences["headers_run2"] == [
            "name",
            "id",
            "value",
        ]
        assert len(order_discrepancy.order_differences["reordered_columns"]) == 2

    def test_column_order_no_difference_same_order(self, tmp_path):
        """Test that no column order discrepancy is reported when order is the same."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        content = "id,name,value\n1,alice,100\n2,bob,200\n"
        file1.write_text(content)
        file2.write_text(content)

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should not find column order difference
        order_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, ColumnOrderDifferenceDiscrepancy)
            ),
            None,
        )
        assert order_discrepancy is None

    def test_column_order_no_difference_different_columns(self, tmp_path):
        """Test that no column order discrepancy is reported when columns are different."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,value\n1,alice,100\n")
        file2.write_text(
            "id,address,phone\n1,123 Main St,555-1234\n"
        )  # Different columns

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should not find column order difference (different columns, not just reordered)
        order_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, ColumnOrderDifferenceDiscrepancy)
            ),
            None,
        )
        assert order_discrepancy is None

        # But should find header difference
        header_discrepancy = next(
            (d for d in discrepancies if isinstance(d, HeaderDifferenceDiscrepancy)),
            None,
        )
        assert header_discrepancy is not None

    def test_no_column_mapping_when_duplicates_exist(self, tmp_path):
        """Test that column mappings are not created when duplicate column names exist."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,id\n1,alice,2\n")  # Duplicate id column
        file2.write_text("id,name,value\n1,alice,100\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        # Should find duplicate column names
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, DuplicateColumnNamesDiscrepancy)
            ),
            None,
        )
        assert dup_discrepancy is not None

        # Should find NO_INDEX_COLUMN discrepancy instead of row differences
        no_index_discrepancy = next(
            (d for d in discrepancies if isinstance(d, NoIndexColumnDiscrepancy)),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy.severity == Severity.CRITICAL

    def test_column_order_difference_detailed_positions(self, tmp_path):
        """Test that column order differences include detailed position information."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("a,b,c,d\n1,2,3,4\n")
        file2.write_text("d,b,a,c\n1,2,3,4\n")  # Completely reordered

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")

        order_discrepancy = next(
            (
                d
                for d in discrepancies
                if isinstance(d, ColumnOrderDifferenceDiscrepancy)
            ),
            None,
        )
        assert order_discrepancy is not None

        reordered_columns = order_discrepancy.order_differences["reordered_columns"]

        # All columns except 'b' should be reordered
        assert len(reordered_columns) == 3  # a, c, d moved positions

        # Check specific position changes
        column_positions = {
            col["column"]: (col["position_run1"], col["position_run2"])
            for col in reordered_columns
        }

        assert column_positions["a"] == (0, 2)  # a moved from position 0 to 2
        assert column_positions["c"] == (2, 3)  # c moved from position 2 to 3
        assert column_positions["d"] == (3, 0)  # d moved from position 3 to 0
        # b stayed at position 1, so it's not in reordered_columns
