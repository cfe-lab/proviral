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

from scripts.compare_runs import (
    Severity,
    Confidence,
    DiscrepancyType,
    Discrepancy,
    ComparisonReport,
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


class TestDiscrepancy:
    """Test the Discrepancy class functionality."""

    def test_discrepancy_creation(self):
        """Test basic discrepancy creation."""
        discrepancy = Discrepancy(
            DiscrepancyType.ROW_DIFFERENCE,
            Severity.HIGH,
            Confidence.HIGH,
            "Test discrepancy",
            location={
                "file": "test.csv",
                "version": "version_7.15",
                "run": "both",
                "row": 2,
            },
            values={"run1": ["A", "B"], "run2": ["A", "C"]},
        )

        assert discrepancy.type == DiscrepancyType.ROW_DIFFERENCE
        assert discrepancy.severity == Severity.HIGH
        assert discrepancy.confidence == Confidence.HIGH
        assert discrepancy.description == "Test discrepancy"
        assert discrepancy.location["file"] == "test.csv"
        assert discrepancy.values["run1"] == ["A", "B"]

    def test_discrepancy_to_dict(self):
        """Test conversion of discrepancy to dictionary."""
        discrepancy = Discrepancy(
            DiscrepancyType.MISSING_FILE,
            Severity.CRITICAL,
            Confidence.HIGH,
            "File missing",
            location={"file": "missing.csv", "version": "version_7.15", "run": "run1"},
            values={"missing_from": "run1"},
        )

        result = discrepancy.to_dict()

        assert result["type"] == "missing_file"
        assert result["severity"] == "CRITICAL"
        assert result["confidence"] == "HIGH"
        assert result["description"] == "File missing"
        assert result["location"]["file"] == "missing.csv"
        assert result["values"]["missing_from"] == "run1"


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
        discrepancy = Discrepancy(
            DiscrepancyType.ROW_DIFFERENCE,
            Severity.MEDIUM,
            Confidence.HIGH,
            "Row differs",
            location={
                "file": "test.csv",
                "version": "version_7.15",
                "run": "both",
                "row": 2,
            },
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

        # Add various discrepancies
        critical_discrepancy = Discrepancy(
            DiscrepancyType.MISSING_FILE, Severity.CRITICAL, Confidence.HIGH, "Critical"
        )
        high_discrepancy = Discrepancy(
            DiscrepancyType.HEADER_DIFFERENCE, Severity.HIGH, Confidence.HIGH, "High"
        )
        medium_discrepancy = Discrepancy(
            DiscrepancyType.ROW_DIFFERENCE, Severity.MEDIUM, Confidence.MEDIUM, "Medium"
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
            (
                d
                for d in discrepancies
                if d.type == DiscrepancyType.ROW_COUNT_DIFFERENCE
            ),
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
            (d for d in discrepancies if d.type == DiscrepancyType.HEADER_DIFFERENCE),
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
            (d for d in discrepancies if d.type == DiscrepancyType.ROW_DIFFERENCE), None
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
            (d for d in discrepancies if d.type == DiscrepancyType.FILE_READ_ERROR),
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
        assert discrepancy["type"] == "missing_file"
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
            "scripts.compare_runs.compare_proviral_files",
            side_effect=Exception("Test error"),
        ):
            report = compare_runs(run1_dir, run2_dir)

        # Should have comparison error discrepancy
        assert "comparison_error" in report.results["version_7.15"]
        discrepancy = report.results["version_7.15"]["comparison_error"][
            "discrepancies"
        ][0]
        assert discrepancy["type"] == "file_read_error"
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
        """Test location fields for row difference discrepancies."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue1\n")
        file2.write_text("header\nvalue2\n")

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")
        row_discrepancy = next(
            (d for d in discrepancies if d.type == DiscrepancyType.ROW_DIFFERENCE), None
        )

        location = row_discrepancy.location

        # Check all required fields are present
        assert "file" in location
        assert "version" in location
        assert "run" in location
        assert "row" in location  # Type-specific field

        assert location["file"] == "test.csv"
        assert location["version"] == "version_7.15"
        assert location["run"] == "both"
        assert location["row"] == 2

    def test_file_read_error_location_fields(self, tmp_path):
        """Test location fields for file read error discrepancies."""
        file1 = tmp_path / "file1.csv"
        file2 = tmp_path / "file2.csv"
        file1.write_text("header\nvalue\n")
        file2.write_text("")  # Empty file

        discrepancies = compare_csv_contents(file1, file2, "version_7.15", "test.csv")
        error_discrepancy = next(
            (d for d in discrepancies if d.type == DiscrepancyType.FILE_READ_ERROR),
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

        with patch("sys.argv", ["compare_runs.py", str(nonexistent), str(existing)]):
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

        with patch("sys.argv", ["compare_runs.py", str(run1_dir), str(run2_dir)]):
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
            ["compare_runs.py", str(run1_dir), str(run2_dir), "-o", str(output_file)],
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
            "sys.argv", ["compare_runs.py", str(run1_dir), str(run2_dir), "--compact"]
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
            "scripts.compare_runs.compare_runs", side_effect=KeyboardInterrupt()
        ):
            with patch("sys.argv", ["compare_runs.py", str(run1_dir), str(run2_dir)]):
                with pytest.raises(SystemExit) as exc_info:
                    main()
                assert exc_info.value.code == 1
