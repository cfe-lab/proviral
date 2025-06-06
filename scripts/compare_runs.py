#! /usr/bin/env python

"""
ws-compare-runs: Compare results from different proviral analysis pipeline runs

This script compares proviral analysis results between different pipeline versions
by analyzing key output files including outcome summaries, proviral landscapes,
table precursors, and primer analysis results.

Examples:
    ws-compare-runs /path/to/run1 /path/to/run2
"""

import argparse
import csv  # must not use pandas or numpy to avoid dependencies.
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
from enum import Enum


# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class Severity(Enum):
    """Severity levels for discrepancies."""

    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


class Confidence(Enum):
    """Confidence levels for discrepancies."""

    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


class DiscrepancyType(Enum):
    """Types of discrepancies that can be detected."""

    MISSING_FILE = "missing_file"
    MISSING_DIRECTORY = "missing_directory"
    HEADER_DIFFERENCE = "header_difference"
    ROW_DIFFERENCE = "row_difference"
    ROW_COUNT_DIFFERENCE = "row_count_difference"
    COLUMN_COUNT_DIFFERENCE = "column_count_difference"
    FILE_READ_ERROR = "file_read_error"
    EMPTY_FILE = "empty_file"


class Discrepancy:
    """Represents a single discrepancy between two runs."""

    def __init__(
        self,
        discrepancy_type: DiscrepancyType,
        severity: Severity,
        confidence: Confidence,
        description: str,
        location: Optional[Dict[str, Any]] = None,
        values: Optional[Dict[str, Any]] = None,
    ):
        self.type = discrepancy_type
        self.severity = severity
        self.confidence = confidence
        self.description = description
        self.location = location or {}
        self.values = values or {}

    def to_dict(self) -> Dict[str, Any]:
        """Convert discrepancy to dictionary for JSON output."""
        return {
            "type": self.type.value,
            "severity": self.severity.value,
            "confidence": self.confidence.value,
            "description": self.description,
            "location": self.location,
            "values": self.values,
        }


class ComparisonReport:
    """Collects and manages all discrepancies found during comparison."""

    def __init__(self, run1_dir: Path, run2_dir: Path):
        self.run1_dir = run1_dir
        self.run2_dir = run2_dir
        self.timestamp = datetime.now().isoformat()
        self.versions_run1: List[str] = []
        self.versions_run2: List[str] = []
        self.common_versions: List[str] = []
        self.results: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(
            lambda: defaultdict(dict)
        )

    def add_discrepancy(self, version: str, filename: str, discrepancy: Discrepancy):
        """Add a discrepancy to the report."""
        if filename not in self.results[version]:
            self.results[version][filename] = {"status": "differs", "discrepancies": []}
        self.results[version][filename]["discrepancies"].append(discrepancy.to_dict())

    def mark_file_identical(self, version: str, filename: str):
        """Mark a file as identical between runs."""
        self.results[version][filename] = {"status": "identical", "discrepancies": []}

    def get_summary(self) -> Dict[str, int]:
        """Generate summary statistics."""
        summary = {
            "total_discrepancies": 0,
            "critical_discrepancies": 0,
            "high_discrepancies": 0,
            "medium_discrepancies": 0,
            "low_discrepancies": 0,
        }

        for version_results in self.results.values():
            for file_results in version_results.values():
                for discrepancy in file_results.get("discrepancies", []):
                    summary["total_discrepancies"] += 1
                    severity = discrepancy["severity"].lower()
                    summary_key = f"{severity}_discrepancies"
                    if summary_key in summary:
                        summary[summary_key] += 1

        return summary

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary for JSON output."""
        return {
            "metadata": {
                "timestamp": self.timestamp,
                "run1_dir": str(self.run1_dir),
                "run2_dir": str(self.run2_dir),
                "versions_run1": self.versions_run1,
                "versions_run2": self.versions_run2,
                "common_versions": self.common_versions,
            },
            "summary": self.get_summary(),
            "results": dict(self.results),
        }

    def to_json(self, indent: int = 2) -> str:
        """Convert report to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)


def find_versions(run_dir: Path) -> List[str]:
    """
    Find all version directories under Results/ subfolder.

    Args:
        run_dir: Path to the run directory

    Returns:
        List of version directory names (e.g., ['version_7.15', 'version_7.17'])
    """
    results_dir = run_dir / "Results"
    if not results_dir.exists():
        logger.warning(f"Results directory not found in {run_dir}")
        return []

    versions = []
    for item in results_dir.iterdir():
        if item.is_dir() and item.name.startswith("version_"):
            versions.append(item.name)

    versions.sort()  # Sort for consistent ordering
    logger.debug(f"Found versions in {run_dir}: {versions}")
    return versions


def get_proviral_csv_files(results_dir: Path, version: str) -> Dict[str, Path]:
    """
    Get all CSV files from the proviral subfolder for a specific version.

    Args:
        results_dir: Path to Results directory
        version: Version directory name

    Returns:
        Dictionary mapping CSV filename to full path
    """
    proviral_dir = results_dir / version / "proviral"
    csv_files: Dict[str, Path] = {}

    if not proviral_dir.exists():
        logger.warning(f"Proviral directory not found: {proviral_dir}")
        return csv_files

    for csv_file in proviral_dir.glob("*.csv"):
        csv_files[csv_file.name] = csv_file

    return csv_files


def read_csv_file(csv_path: Path) -> List[List[str]]:
    """
    Read a CSV file and return its contents as a list of lists.

    Args:
        csv_path: Path to the CSV file

    Returns:
        List of rows, where each row is a list of values
    """
    try:
        with open(csv_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            return list(reader)
    except Exception as e:
        logger.error(f"Error reading CSV file {csv_path}: {e}")
        return []


def compare_csv_contents(
    file1: Path, file2: Path, version: str, filename: str
) -> List[Discrepancy]:
    """
    Compare the contents of two CSV files and return list of discrepancies.

    Args:
        file1: Path to first CSV file
        file2: Path to second CSV file
        version: Version being compared
        filename: Name of the file being compared

    Returns:
        List of Discrepancy objects
    """
    content1 = read_csv_file(file1)
    content2 = read_csv_file(file2)

    discrepancies: List[Discrepancy] = []

    # Check if either file is empty or failed to read
    if not content1 and not content2:
        return discrepancies  # Both empty, no discrepancies

    if not content1:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.FILE_READ_ERROR,
                Severity.CRITICAL,
                Confidence.HIGH,
                "First file is empty or failed to read",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run1",
                    "file_path": str(file1),
                },
            )
        )
        return discrepancies

    if not content2:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.FILE_READ_ERROR,
                Severity.CRITICAL,
                Confidence.HIGH,
                "Second file is empty or failed to read",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run2",
                    "file_path": str(file2),
                },
            )
        )
        return discrepancies

    # Check if files have same number of rows
    if len(content1) != len(content2):
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.ROW_COUNT_DIFFERENCE,
                Severity.HIGH,
                Confidence.HIGH,
                f"Row count differs: {len(content1)} vs {len(content2)}",
                location={"file": filename, "version": version, "run": "both"},
                values={"run1_rows": len(content1), "run2_rows": len(content2)},
            )
        )

    # Compare row by row
    max_rows = max(len(content1), len(content2))
    for i in range(max_rows):
        row1 = content1[i] if i < len(content1) else []
        row2 = content2[i] if i < len(content2) else []

        if row1 != row2:
            if i == 0:  # Header row
                discrepancies.append(
                    Discrepancy(
                        DiscrepancyType.HEADER_DIFFERENCE,
                        Severity.CRITICAL,
                        Confidence.HIGH,
                        "Header row differs",
                        location={
                            "file": filename,
                            "version": version,
                            "run": "both",
                            "row": i + 1,
                        },
                        values={"run1": row1, "run2": row2},
                    )
                )

                # Check for column count differences in header
                if len(row1) != len(row2):
                    discrepancies.append(
                        Discrepancy(
                            DiscrepancyType.COLUMN_COUNT_DIFFERENCE,
                            Severity.CRITICAL,
                            Confidence.HIGH,
                            f"Header column count differs: {len(row1)} vs {len(row2)}",
                            location={
                                "file": filename,
                                "version": version,
                                "run": "both",
                                "row": i + 1,
                            },
                            values={
                                "run1_columns": len(row1),
                                "run2_columns": len(row2),
                            },
                        )
                    )
            else:
                # Determine severity based on content analysis
                severity = _determine_row_difference_severity(row1, row2, i)
                confidence = _determine_row_difference_confidence(row1, row2)

                discrepancies.append(
                    Discrepancy(
                        DiscrepancyType.ROW_DIFFERENCE,
                        severity,
                        confidence,
                        f"Row {i + 1} differs",
                        location={
                            "file": filename,
                            "version": version,
                            "run": "both",
                            "row": i + 1,
                        },
                        values={"run1": row1, "run2": row2},
                    )
                )

                # Check for column count differences in data rows
                if len(row1) != len(row2):
                    discrepancies.append(
                        Discrepancy(
                            DiscrepancyType.COLUMN_COUNT_DIFFERENCE,
                            Severity.HIGH,
                            Confidence.HIGH,
                            f"Row {i + 1} column count differs: {len(row1)} vs {len(row2)}",
                            location={
                                "file": filename,
                                "version": version,
                                "run": "both",
                                "row": i + 1,
                            },
                            values={
                                "run1_columns": len(row1),
                                "run2_columns": len(row2),
                            },
                        )
                    )

    return discrepancies


def _determine_row_difference_severity(
    row1: List[str], row2: List[str], row_index: int
) -> Severity:
    """Determine severity of row difference based on content analysis."""
    # For proviral analysis, certain columns are more critical
    critical_outcomes = {"failed", "success", "complete", "incomplete"}

    # Check if key outcome fields differ
    for val1, val2 in zip(row1, row2):
        if val1.lower() in critical_outcomes or val2.lower() in critical_outcomes:
            if val1 != val2:
                return Severity.HIGH

    # Check if any numeric values differ significantly
    for val1, val2 in zip(row1, row2):
        try:
            num1, num2 = float(val1), float(val2)
            if abs(num1 - num2) > 0.001:  # Significant numeric difference
                return Severity.MEDIUM
        except ValueError:
            continue

    # Default to medium for any data difference
    return Severity.MEDIUM


def _determine_row_difference_confidence(
    row1: List[str], row2: List[str]
) -> Confidence:
    """Determine confidence level for row difference."""
    # High confidence if rows are clearly different
    if len(row1) != len(row2):
        return Confidence.HIGH

    # Check for clear value differences
    different_values = sum(1 for v1, v2 in zip(row1, row2) if v1 != v2)
    if different_values > len(row1) / 2:  # More than half the values differ
        return Confidence.HIGH
    elif different_values > 1:  # Multiple values differ
        return Confidence.HIGH
    else:  # Single value difference
        return Confidence.MEDIUM


def compare_proviral_files(
    run1_dir: Path, run2_dir: Path, version: str, report: ComparisonReport
) -> None:
    """
    Compare proviral CSV files between two runs for a specific version.

    Args:
        run1_dir: Path to first run directory
        run2_dir: Path to second run directory
        version: Version to compare
        report: ComparisonReport to add discrepancies to
    """
    logger.debug(f"Comparing proviral files for {version}")

    # Get CSV files from both runs
    csv_files1 = get_proviral_csv_files(run1_dir / "Results", version)
    csv_files2 = get_proviral_csv_files(run2_dir / "Results", version)

    # Find all unique CSV file names
    all_csv_names = set(csv_files1.keys()) | set(csv_files2.keys())

    if not all_csv_names:
        report.add_discrepancy(
            version,
            "proviral_directory",
            Discrepancy(
                DiscrepancyType.MISSING_DIRECTORY,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"No proviral CSV files found for {version}",
                location={
                    "file": "proviral_directory",
                    "version": version,
                    "run": "both",
                },
            ),
        )
        return

    for csv_name in sorted(all_csv_names):
        if csv_name not in csv_files1:
            report.add_discrepancy(
                version,
                csv_name,
                Discrepancy(
                    DiscrepancyType.MISSING_FILE,
                    Severity.CRITICAL,
                    Confidence.HIGH,
                    f"File missing in run1: {csv_name}",
                    location={
                        "file": csv_name,
                        "version": version,
                        "run": "run1",
                        "missing_from": "run1",
                    },
                ),
            )
            continue

        if csv_name not in csv_files2:
            report.add_discrepancy(
                version,
                csv_name,
                Discrepancy(
                    DiscrepancyType.MISSING_FILE,
                    Severity.CRITICAL,
                    Confidence.HIGH,
                    f"File missing in run2: {csv_name}",
                    location={
                        "file": csv_name,
                        "version": version,
                        "run": "run2",
                        "missing_from": "run2",
                    },
                ),
            )
            continue

        # Compare the files
        discrepancies = compare_csv_contents(
            csv_files1[csv_name], csv_files2[csv_name], version, csv_name
        )

        if not discrepancies:
            report.mark_file_identical(version, csv_name)
        else:
            for discrepancy in discrepancies:
                report.add_discrepancy(version, csv_name, discrepancy)


def compare_runs(run1_dir: Path, run2_dir: Path) -> ComparisonReport:
    """
    Compare proviral analysis results between two runs.

    Args:
        run1_dir: Path to first run directory
        run2_dir: Path to second run directory

    Returns:
        ComparisonReport with all discrepancies found
    """
    logger.debug(f"Comparing runs: {run1_dir} vs {run2_dir}")

    report = ComparisonReport(run1_dir, run2_dir)

    # Find versions in both runs
    versions1 = find_versions(run1_dir)
    versions2 = find_versions(run2_dir)

    report.versions_run1 = versions1
    report.versions_run2 = versions2

    # Find common versions
    common_versions = set(versions1) & set(versions2)
    report.common_versions = sorted(common_versions)

    if not common_versions:
        # Add discrepancy for no common versions
        report.add_discrepancy(
            "metadata",
            "versions",
            Discrepancy(
                DiscrepancyType.MISSING_DIRECTORY,
                Severity.CRITICAL,
                Confidence.HIGH,
                "No common versions found between the two runs",
                location={"file": "versions", "version": "all", "run": "both"},
                values={"run1_versions": versions1, "run2_versions": versions2},
            ),
        )
        return report

    # Report version mismatches
    only_in_run1 = set(versions1) - set(versions2)
    only_in_run2 = set(versions2) - set(versions1)

    if only_in_run1:
        report.add_discrepancy(
            "metadata",
            "versions",
            Discrepancy(
                DiscrepancyType.MISSING_DIRECTORY,
                Severity.MEDIUM,
                Confidence.HIGH,
                f"Versions only in run1: {sorted(only_in_run1)}",
                location={"file": "versions", "version": "multiple", "run": "run1"},
                values={"missing_versions": sorted(only_in_run1)},
            ),
        )

    if only_in_run2:
        report.add_discrepancy(
            "metadata",
            "versions",
            Discrepancy(
                DiscrepancyType.MISSING_DIRECTORY,
                Severity.MEDIUM,
                Confidence.HIGH,
                f"Versions only in run2: {sorted(only_in_run2)}",
                location={"file": "versions", "version": "multiple", "run": "run2"},
                values={"missing_versions": sorted(only_in_run2)},
            ),
        )

    # Compare each common version
    for version in sorted(common_versions):
        try:
            compare_proviral_files(run1_dir, run2_dir, version, report)
        except Exception as e:
            report.add_discrepancy(
                version,
                "comparison_error",
                Discrepancy(
                    DiscrepancyType.FILE_READ_ERROR,
                    Severity.CRITICAL,
                    Confidence.MEDIUM,
                    f"Error comparing {version}: {str(e)}",
                    location={
                        "file": "comparison_error",
                        "version": version,
                        "run": "both",
                        "error": str(e),
                    },
                ),
            )

    return report


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Compare proviral analysis results between different pipeline runs",
        epilog="Example: ws-compare-runs /path/to/run1 /path/to/run2",
    )

    parser.add_argument("run1_dir", type=Path, help="Path to first run directory")

    parser.add_argument("run2_dir", type=Path, help="Path to second run directory")

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    parser.add_argument(
        "-o", "--output", type=Path, help="Output JSON file (default: stdout)"
    )

    parser.add_argument(
        "--compact", action="store_true", help="Output compact JSON (no indentation)"
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        # Reduce logging to only show errors
        logging.getLogger().setLevel(logging.ERROR)

    # Validate input directories
    if not args.run1_dir.exists():
        logger.error(f"Run1 directory does not exist: {args.run1_dir}")
        sys.exit(1)

    if not args.run2_dir.exists():
        logger.error(f"Run2 directory does not exist: {args.run2_dir}")
        sys.exit(1)

    if not args.run1_dir.is_dir():
        logger.error(f"Run1 path is not a directory: {args.run1_dir}")
        sys.exit(1)

    if not args.run2_dir.is_dir():
        logger.error(f"Run2 path is not a directory: {args.run2_dir}")
        sys.exit(1)

    # Perform the comparison
    try:
        report = compare_runs(args.run1_dir, args.run2_dir)

        # Generate JSON output
        indent = None if args.compact else 2
        json_output = report.to_json(indent=indent)

        # Output JSON
        if args.output:
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(json_output)
            if args.verbose:
                print(f"Report written to {args.output}")
        else:
            print(json_output)

    except KeyboardInterrupt:
        logger.error("Comparison interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error during comparison: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
