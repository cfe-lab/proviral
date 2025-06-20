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
    compare_runs,
    compare_proviral_files,
    compare_csv_contents,
    main,
    _determine_row_difference_severity,
    _determine_row_difference_confidence,
)
from scripts.compare_runs.file_operations import (
    find_versions,
    get_proviral_csv_files,
    read_csv_file,
)
from scripts.compare_runs.severity import Severity
from scripts.compare_runs.confidence import Confidence
from scripts.compare_runs.comparison_report import ComparisonReport

# Import all discrepancy class types and error types referenced in tests
from scripts.compare_runs.discrepancy import (
    MissingFile,
    MissingDirectory,
    HeaderDifference,
    RowDifference,
    RowCountDifference,
    ColumnCountDifference,
    DuplicateColumnNames,
    ColumnOrderDifference,
    RowOrderDifference,
    MissingRow,
    ExtraRow,
    FieldChange,
    HeaderFieldChange,
    ColumnReorder,
    RowReorder,
)
from scripts.compare_runs.errors import (
    FileReadError,
    NoIndexColumn,
)


# Helper function to wrap the new compare_csv_contents for backward compatibility
def legacy_compare_csv_contents(
    file1: Path,
    file2: Path,
    version: str,
    filename: str,
    index_pattern: str = ".*/header1",
):
    """Wrapper for compare_csv_contents that returns discrepancies like the old signature."""
    report = ComparisonReport(file1.parent, file2.parent)
    compare_csv_contents(report, file1, file2, version, filename, index_pattern)
    # Return both discrepancies and errors as one list for backward compatibility
    all_items = []
    all_items.extend([discrepancy.to_dict() for discrepancy in report.results])
    all_items.extend([error.to_dict() for error in report.errors])
    return all_items


class TestDiscrepancy:
    """Test the Discrepancy class functionality."""

    def test_discrepancy_creation(self):
        """Test basic discrepancy creation."""
        # Changed to use a specific discrepancy type: RowDifferenceDiscrepancy
        discrepancy = RowDifference(
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

        assert isinstance(discrepancy, RowDifference)
        assert discrepancy.severity == Severity.HIGH
        assert discrepancy.confidence == Confidence.HIGH
        assert discrepancy.description == "Test discrepancy"
        # Assertions for location and values properties (should still work due to base class implementation)
        result_dict = discrepancy.to_dict()
        assert result_dict["location"]["file"] == "test.csv"
        assert result_dict["location"]["row"] == 2

    def test_discrepancy_to_dict(self):
        """Test conversion of discrepancy to dictionary."""
        # Changed to use a specific discrepancy type: MissingFileDiscrepancy
        discrepancy = MissingFile(
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
        assert result["location"]["missing_from"] == "run1"


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
        discrepancy = RowDifference(
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

        report.add_discrepancy(discrepancy)

        assert len([d.to_dict() for d in report.results]) == 1
        assert [d.to_dict() for d in report.results][0][
            "type"
        ] == "RowDifference"  # Trimmed by _trim_data_recursively
        assert [d.to_dict() for d in report.results][0]["location"][
            "version"
        ] == "version_7.15"
        assert [d.to_dict() for d in report.results][0]["location"][
            "file"
        ] == "test.csv"

    def test_mark_file_identical(self):
        """Test marking files as identical."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))

        report.mark_file_identical()

        # In the flat structure, identical files don't generate discrepancies
        assert len([d.to_dict() for d in report.results]) == 0

    def test_get_summary(self):
        """Test summary statistics generation."""
        report = ComparisonReport(Path("/fake/run1"), Path("/fake/run2"))

        # Add various discrepancies using specific types
        critical_discrepancy = MissingFile(
            severity=Severity.CRITICAL,
            confidence=Confidence.HIGH,
            description="Critical",
            file="file1.csv",
            version="version_7.15",
            run="run1",  # run where it's missing
            missing_from="run1",
            present_in="run2",  # Placeholders
        )
        high_discrepancy = HeaderDifference(
            severity=Severity.HIGH,
            confidence=Confidence.HIGH,
            description="High",
            file="file2.csv",
            version="version_7.15",
            run="both",  # Comparison of both
            # Placeholders for HeaderDifferenceDiscrepancy fields
            columns_missing_in_run2=[],
            columns_extra_in_run2=[],
            run1_header_count=3,
            run2_header_count=3,
        )
        medium_discrepancy = RowDifference(
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

        report.add_discrepancy(critical_discrepancy)
        report.add_discrepancy(high_discrepancy)
        report.add_discrepancy(medium_discrepancy)

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
        report = ComparisonReport(tmp_path / "run1", tmp_path / "run2")

        compare_csv_contents(report, file1, file2, "version_7.15", "test.csv", ".*/.*")

        assert [d.to_dict() for d in report.results] == []

    def test_compare_different_row_count(self, tmp_path):
        """Test comparison with different row counts."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header2\nvalue1,value2\nvalue3,value4\n")
        report = ComparisonReport(tmp_path / "run1", tmp_path / "run2")

        compare_csv_contents(
            report, file1, file2, "version_7.15", "test.csv", ".*/header1"
        )
        discrepancies = [d.to_dict() for d in report.results]

        # Should find row count difference
        row_count_discrepancy = next(
            (
                d
                for d in discrepancies
                if d.get("type") == "RowCountDifference"
            ),
            None,
        )
        assert row_count_discrepancy is not None
        assert row_count_discrepancy["severity"] == "HIGH"
        result_dict = row_count_discrepancy
        assert result_dict["location"]["file"] == "test.csv"
        assert result_dict["location"]["version"] == "version_7.15"
        assert result_dict["location"]["run"] == "both"

    def test_compare_different_headers(self, tmp_path):
        """Test comparison with different headers."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header3\nvalue1,value2\n")

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/header1"
        )

        # Should find header field changes
        header_discrepancies = [
            d for d in discrepancies if d.get("type") == "HeaderFieldChange"
        ]
        assert len(header_discrepancies) == 1  # header2 -> header3 change

        header_discrepancy = header_discrepancies[0]
        assert header_discrepancy["severity"] == "CRITICAL"
        result_dict = header_discrepancy
        assert result_dict["location"]["file"] == "test.csv"
        assert result_dict["location"]["row"] == 1
        assert result_dict["location"]["column_index"] == 1  # Second column
        assert result_dict["value1"] == "header2"
        assert result_dict["value2"] == "header3"

    def test_compare_different_data_rows(self, tmp_path):
        """Test comparison with different data rows."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\nvalue1,value2\n")
        file2.write_text("header1,header2\nvalue1,different\n")

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find field change discrepancies
        field_discrepancies = [
            d for d in discrepancies if d.get("type") == "FieldChange"
        ]
        assert len(field_discrepancies) == 1  # One field changed: header2 column

        field_discrepancy = field_discrepancies[0]
        result_dict = field_discrepancy
        assert result_dict["location"]["file"] == "test.csv"
        assert result_dict["location"]["row"] == 2
        assert result_dict["location"]["column_index"] == 1  # Second column
        assert result_dict["value1"] == "value2"
        assert result_dict["value2"] == "different"

    def test_compare_empty_files(self, tmp_path):
        """Test comparison of empty files."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("")
        file2.write_text("")

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        assert discrepancies == []

    def test_compare_one_empty_file(self, tmp_path):
        """Test comparison where one file is empty."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2\n")
        file2.write_text("")

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find file read error
        error_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "FileReadError"),
            None,
        )
        assert error_discrepancy is not None
        assert error_discrepancy["severity"] == "CRITICAL"
        result_dict = error_discrepancy
        assert result_dict["location"]["run"] == "run2"


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

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report, ".*/.*")

        # Should have missing directory discrepancy
        assert len([d.to_dict() for d in report.results]) > 0
        assert any(
            d["location"]["file"] == "proviral_directory"
            and d["location"]["version"] == "version_7.15"
            for d in [d.to_dict() for d in report.results]
        )

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

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report, ".*/.*")

        # Should have missing file discrepancy
        assert len([d.to_dict() for d in report.results]) > 0
        missing_file_discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "missing.csv"
            ),
            None,
        )
        assert missing_file_discrepancy is not None
        assert missing_file_discrepancy["location"]["missing_from"] == "run1"

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

        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report, ".*/.*")

        # Should have no discrepancies for identical files
        identical_file_discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "identical.csv"
            ),
            None,
        )
        assert identical_file_discrepancy is None


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

        report = compare_runs(run1_dir, run2_dir, ".*/.*")

        # Should have no common versions discrepancy
        assert len([d.to_dict() for d in report.results]) > 0
        no_common_versions_discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "versions"
                and d["location"]["version"] == "all"
            ),
            None,
        )
        assert no_common_versions_discrepancy is not None

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

        report = compare_runs(run1_dir, run2_dir, ".*/.*")

        assert report.common_versions == ["version_7.15"]
        # Should have discrepancies for version mismatches
        version_mismatch_discrepancies = [
            d
            for d in [d.to_dict() for d in report.results]
            if d["location"]["file"] == "versions"
        ]
        assert len(version_mismatch_discrepancies) > 0

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

        report = compare_runs(run1_dir, run2_dir, ".*/.*")

        assert report.common_versions == ["version_7.15"]
        # Should have errors for the files that can't be compared due to no index column
        # or discrepancies in results list for actual differences found
        has_results_or_errors = len(report.results) > 0 or len(report.errors) > 0
        assert has_results_or_errors

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
            report = compare_runs(run1_dir, run2_dir, ".*/.*")

        # Should have comparison error discrepancy in errors list
        comparison_error_discrepancy = next(
            (d for d in report.errors if "Test error" in d.description),
            None,
        )
        assert comparison_error_discrepancy is not None
        assert "Test error" in comparison_error_discrepancy.description


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
        compare_proviral_files(run1_dir, run2_dir, "version_7.15", report, ".*/.*")

        discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "missing.csv"
            ),
            None,
        )
        assert discrepancy is not None
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
        """Test that NO_INDEX_COLUMN discrepancies are reported when regex doesn't match any columns."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue1\n")
        file2.write_text("header\nvalue2\n")

        # Use a regex that doesn't match any columns to force NO_INDEX_COLUMN error
        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", "nomatch/.*"
        )

        # Should find NO_INDEX_COLUMN discrepancy instead of row differences
        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
            None,
        )
        assert no_index_discrepancy is not None

        result_dict = no_index_discrepancy
        location = result_dict["location"]

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "reason" in result_dict

        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"

    def test_file_read_error_location_fields(self, tmp_path):
        """Test location fields for file read error discrepancies."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue\n")
        file2.write_text("")  # Empty file

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )
        error_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "FileReadError"),
            None,
        )

        result_dict = error_discrepancy
        location = result_dict["location"]

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "file_path" in location

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

        with patch(
            "sys.argv",
            ["compare_runs", "--indexes", ".*/.*", str(nonexistent), str(existing)],
        ):
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

        with patch(
            "sys.argv",
            ["compare_runs", "--indexes", ".*/.*", str(run1_dir), str(run2_dir)],
        ):
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
            [
                "compare_runs",
                "--indexes",
                ".*/.*",
                str(run1_dir),
                str(run2_dir),
                "-o",
                str(output_file),
            ],
        ):
            main()

        assert output_file.exists()
        output = json.loads(output_file.read_text())
        assert "metadata" in output

    def test_main_with_comact_output(self, tmp_path, capsys):
        """Test main function with compact JSON output."""
        # Setup directories
        run1_dir = tmp_path / "run1"
        run2_dir = tmp_path / "run2"

        for run_dir in [run1_dir, run2_dir]:
            proviral_dir = run_dir / "Results" / "version_7.15" / "proviral"
            proviral_dir.mkdir(parents=True)
            (proviral_dir / "test.csv").write_text("header\nvalue\n")

        with patch(
            "sys.argv",
            [
                "compare_runs",
                "--indexes",
                ".*/.*",
                str(run1_dir),
                str(run2_dir),
                "--compact",
            ],
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
            with patch(
                "sys.argv",
                ["compare_runs", "--indexes", ".*/.*", str(run1_dir), str(run2_dir)],
            ):
                with pytest.raises(SystemExit) as exc_info:
                    main()
                assert exc_info.value.code == 1


class TestCSVComparisonWithIndexDiscovery:
    """Test CSV comparison functionality with index column discovery."""

    def test_compare_csv_with_index_discovery_identical_files(self, tmp_path):
        """Test CSV comparison with index discovery for identical files."""
        content = "sample_id,result\nS001,pass\nS002,fail\n"
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text(content)
        file2.write_text(content)

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/sample_id"
        )

        # Should be no discrepancies for identical files
        assert discrepancies == []

    def test_compare_csv_with_no_index_column_discrepancy(self, tmp_path):
        """Test CSV comparison generates NO_INDEX_COLUMN discrepancy when no columns match the regex pattern."""
        # Create files with columns that won't match a specific pattern
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,value\n1,A\n2,B\n")
        file2.write_text("id,value\n3,C\n4,D\n")

        # Use a regex that won't match any columns
        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", "nonexistent/.*"
        )

        # Should find NO_INDEX_COLUMN discrepancy
        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy["severity"] == "CRITICAL"
        assert no_index_discrepancy["confidence"] == "HIGH"
        assert "Row comparison skipped" in no_index_discrepancy["description"]

    def test_compare_csv_no_index_discrepancy_when_no_unique_columns(self, tmp_path):
        """Test that NO_INDEX_COLUMN discrepancy is generated when regex doesn't match any columns."""
        # Create files with duplicate columns - doesn't matter for regex matching
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text(
            "category,status\nA,active\nA,active\n"
        )  # Both columns have duplicates
        file2.write_text(
            "category,status\nB,pending\nB,pending\n"
        )  # Both columns have duplicates

        # Use a regex that doesn't match any columns
        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", "missing_file/.*"
        )

        # Should find NO_INDEX_COLUMN discrepancy since no columns match the pattern
        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy["severity"] == "CRITICAL"
        assert "Row comparison skipped" in no_index_discrepancy["description"]

    def test_compare_csv_with_successful_index_discovery(self, tmp_path):
        """Test CSV comparison with successful index column discovery."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("sample_id,result\nS001,pass\nS002,fail\n")
        file2.write_text(
            "sample_id,result\nS001,pass\nS002,pass\n"  # S002 result differs
        )

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/sample_id"
        )

        # Should find field change discrepancy but no NO_INDEX_COLUMN discrepancy
        field_discrepancies = [
            d for d in discrepancies if d.get("type") == "FieldChange"
        ]
        assert len(field_discrepancies) == 1  # One field changed: S002 result

        field_discrepancy = field_discrepancies[0]
        assert field_discrepancy["value1"] == "fail"
        assert field_discrepancy["value2"] == "pass"

        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
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

        # Use a regex that doesn't match any columns to force NO_INDEX_COLUMN error
        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", "missing/.*"
        )
        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
            None,
        )

        assert no_index_discrepancy is not None
        result_dict = no_index_discrepancy
        location = result_dict["location"]

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
        compare_proviral_files(
            run1_dir, run2_dir, "version_7.15", report, ".*/sample_id"
        )

        # In the flat structure, index info is embedded in discrepancies
        # For files with differences, check if there are discrepancies
        test_file_discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "test.csv"
            ),
            None,
        )
        assert test_file_discrepancy is not None

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
        compare_proviral_files(
            run1_dir, run2_dir, "version_7.15", report, ".*/sample_id"
        )

        # In the flat structure, identical files don't generate discrepancies
        identical_file_discrepancy = next(
            (
                d
                for d in [d.to_dict() for d in report.results]
                if d["location"]["file"] == "identical.csv"
            ),
            None,
        )
        assert identical_file_discrepancy is None


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

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find duplicate column names discrepancy
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if d.get("type") == "DuplicateColumnNames"
            ),
            None,
        )
        assert dup_discrepancy is not None
        assert dup_discrepancy["severity"] == "CRITICAL"
        assert dup_discrepancy["confidence"] == "HIGH"
        assert "Duplicate column name 'header1'" in dup_discrepancy["description"]
        assert "run1" in dup_discrepancy["description"]

        # Check location fields
        result_dict = dup_discrepancy
        location = result_dict["location"]
        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "run1"
        assert result_dict["duplicate_column"] == "header1"
        assert location["positions"] == [0, 2]  # header1 appears at positions 0 and 2

        # Check values
        assert dup_discrepancy["duplicate_column"] == "header1"
        assert dup_discrepancy["all_headers"] == ["header1", "header2", "header1"]

    def test_duplicate_column_names_run2(self, tmp_path):
        """Test detection of duplicate column names in run2."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header1,header2,header3\nvalue1,value2,value3\n")
        file2.write_text(
            "col1,col2,col1,col3\nvalue1,value2,value3,value4\n"
        )  # Duplicate col1

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find duplicate column names discrepancy
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if d.get("type") == "DuplicateColumnNames"
            ),
            None,
        )
        assert dup_discrepancy is not None
        assert "Duplicate column name 'col1'" in dup_discrepancy["description"]
        assert "run2" in dup_discrepancy["description"]
        result_dict = dup_discrepancy
        assert result_dict["location"]["run"] == "run2"
        assert result_dict["duplicate_column"] == "col1"
        assert result_dict["location"]["positions"] == [0, 2]

    def test_duplicate_column_names_both_runs(self, tmp_path):
        """Test detection of duplicate column names in both runs."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,id\nvalue1,alice,value3\n")  # Duplicate id
        file2.write_text("name,id,name\nvalue1,value2,value3\n")  # Duplicate name

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find duplicate column names discrepancies for both runs
        dup_discrepancies = [
            d
            for d in discrepancies
            if d.get("type") == "DuplicateColumnNames"
        ]
        assert len(dup_discrepancies) == 2

        # Check that we have one for each run
        run1_discrepancy = next(
            (d for d in dup_discrepancies if d["location"]["run"] == "run1"),
            None,
        )
        run2_discrepancy = next(
            (d for d in dup_discrepancies if d["location"]["run"] == "run2"),
            None,
        )

        assert run1_discrepancy is not None
        assert run2_discrepancy is not None
        assert run1_discrepancy["duplicate_column"] == "id"
        assert run2_discrepancy["duplicate_column"] == "name"

    def test_column_order_difference(self, tmp_path):
        """Test detection of column order differences."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,value\n1,alice,100\n2,bob,200\n")
        file2.write_text(
            "name,id,value\n1,alice,100\n2,bob,200\n"
        )  # Same columns, different order

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should find column reorder discrepancies
        reorder_discrepancies = [
            d for d in discrepancies if d.get("type") == "ColumnReorder"
        ]
        assert len(reorder_discrepancies) == 2  # id and name both moved

        # Check that we got the expected column reorders
        reorder_cols = {d["location"]["column_name"]: d for d in reorder_discrepancies}
        assert "id" in reorder_cols
        assert "name" in reorder_cols

        # id moved from position 0 to 1
        id_discrepancy = reorder_cols["id"]
        assert id_discrepancy["severity"] == "MEDIUM"
        assert id_discrepancy["confidence"] == "HIGH"
        assert id_discrepancy["location"]["position_run1"] == 0
        assert id_discrepancy["location"]["position_run2"] == 1

        # name moved from position 1 to 0
        name_discrepancy = reorder_cols["name"]
        assert name_discrepancy["location"]["position_run1"] == 1
        assert name_discrepancy["location"]["position_run2"] == 0

        # Check location fields
        result_dict = id_discrepancy
        location = result_dict["location"]
        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"
        assert location["column_name"] == "id"

    def test_column_order_no_difference_same_order(self, tmp_path):
        """Test that no column order discrepancy is reported when order is the same."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        content = "id,name,value\n1,alice,100\n2,bob,200\n"
        file1.write_text(content)
        file2.write_text(content)

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv"
        )

        # Should not find column order difference
        order_discrepancy = next(
            (
                d
                for d in discrepancies
                if d.get("type") == "ColumnOrderDifference"
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

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/id"
        )

        # Should not find column reorder discrepancies (different columns, not just reordered)
        reorder_discrepancies = [
            d for d in discrepancies if d.get("type") == "ColumnReorder"
        ]
        assert len(reorder_discrepancies) == 0

        # But should find header field changes
        header_discrepancies = [
            d for d in discrepancies if d.get("type") == "HeaderFieldChange"
        ]
        assert len(header_discrepancies) >= 1  # At least one header change

    def test_no_column_mapping_when_duplicates_exist(self, tmp_path):
        """Test that column mappings are not created when duplicate column names exist."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("id,name,id\n1,alice,2\n")  # Duplicate id column
        file2.write_text("id,name,value\n1,alice,100\n")

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/name"
        )

        # Should find duplicate column names
        dup_discrepancy = next(
            (
                d
                for d in discrepancies
                if d.get("type") == "DuplicateColumnNames"
            ),
            None,
        )
        assert dup_discrepancy is not None

        # Should find NO_INDEX_COLUMN discrepancy instead of row differences
        no_index_discrepancy = next(
            (d for d in discrepancies if d.get("type") == "NoIndexColumn"),
            None,
        )
        assert no_index_discrepancy is not None
        assert no_index_discrepancy["severity"] == "CRITICAL"

    def test_column_order_difference_detailed_positions(self, tmp_path):
        """Test that column order differences include detailed position information."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("a,b,c,d\n1,2,3,4\n")
        file2.write_text("d,b,a,c\n1,2,3,4\n")  # Completely reordered

        discrepancies = legacy_compare_csv_contents(
            file1, file2, "version_7.15", "test.csv", ".*/a"
        )

        reorder_discrepancies = [
            d for d in discrepancies if d.get("type") == "ColumnReorder"
        ]

        # All columns except 'b' should be reordered
        assert len(reorder_discrepancies) == 3  # a, c, d moved positions

        # Check specific position changes
        column_positions = {
            d["location"]["column_name"]: (
                d["location"]["position_run1"],
                d["location"]["position_run2"],
            )
            for d in reorder_discrepancies
        }

        assert column_positions["a"] == (0, 2)  # a moved from position 0 to 2
        assert column_positions["c"] == (2, 3)  # c moved from position 2 to 3
        assert column_positions["d"] == (3, 0)  # d moved from position 3 to 0
        # b stayed at position 1, so it's not in reordered_columns
