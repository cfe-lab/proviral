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
    NO_INDEX_COLUMN = "no_index_column"
    DUPLICATE_COLUMN_NAMES = "duplicate_column_names"
    COLUMN_ORDER_DIFFERENCE = "column_order_difference"
    ROW_ORDER_DIFFERENCE = "row_order_difference"
    MISSING_ROW = "missing_row"
    EXTRA_ROW = "extra_row"


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
            self.results[version][filename] = {
                "status": "differs",
                "discrepancies": [],
                "index_column": None,
            }
        self.results[version][filename]["discrepancies"].append(discrepancy.to_dict())

    def mark_file_identical(
        self, version: str, filename: str, index_info: Optional[Dict[str, Any]] = None
    ):
        """Mark a file as identical between runs."""
        self.results[version][filename] = {
            "status": "identical",
            "discrepancies": [],
            "index_column": index_info,
        }

    def set_file_index_info(
        self, version: str, filename: str, index_info: Optional[Dict[str, Any]]
    ):
        """Set index column information for a file."""
        if filename not in self.results[version]:
            self.results[version][filename] = {
                "status": "unknown",
                "discrepancies": [],
                "index_column": None,
            }
        self.results[version][filename]["index_column"] = index_info

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


def _analyze_header_differences(
    header1: List[str], header2: List[str]
) -> Dict[str, Any]:
    """Analyze differences between two header rows."""
    changes: Dict[str, Any] = {
        "indices": [],
        "added": [],
        "removed": [],
        "modified": [],
    }

    max_len = max(len(header1), len(header2))

    for i in range(max_len):
        val1 = header1[i] if i < len(header1) else None
        val2 = header2[i] if i < len(header2) else None

        if val1 != val2:
            changes["indices"].append(i)

            if val1 is None:
                changes["added"].append({"index": i, "value": val2})
            elif val2 is None:
                changes["removed"].append({"index": i, "value": val1})
            else:
                changes["modified"].append(
                    {"index": i, "old_header": val1, "new_header": val2}
                )

    return changes


def _get_header_change_summary(changes: Dict[str, Any]) -> str:
    """Generate a human-readable summary of header changes."""
    parts = []

    if changes["added"]:
        parts.append(f"{len(changes['added'])} headers added")
    if changes["removed"]:
        parts.append(f"{len(changes['removed'])} headers removed")
    if changes["modified"]:
        parts.append(f"{len(changes['modified'])} headers modified")

    if not parts:
        return "unknown header changes"

    return ", ".join(parts)


def _analyze_row_differences(
    row1: List[str],
    row2: List[str],
    headers1: List[str],
    headers2: List[str],
    column_map1: Optional[Dict[str, int]] = None,
    column_map2: Optional[Dict[str, int]] = None,
) -> Dict[str, Any]:
    """
    Analyze differences between two data rows using column names when possible.

    Args:
        row1: First row data
        row2: Second row data
        headers1: Headers from first file
        headers2: Headers from second file
        column_map1: Column name to index mapping for first file
        column_map2: Column name to index mapping for second file

    Returns:
        Dictionary with change information including column names
    """
    changes: Dict[str, Any] = {
        "indices": [],
        "change_types": [],
        "field_changes": [],
        "column_names": [],
    }

    # Use column mappings if available, otherwise fall back to index-based comparison
    if (
        column_map1
        and column_map2
        and set(column_map1.keys()) == set(column_map2.keys())
    ):
        # Compare by column name when both files have same column structure
        for col_name in column_map1.keys():
            idx1 = column_map1[col_name]
            idx2 = column_map2[col_name]

            val1 = row1[idx1] if idx1 < len(row1) else ""
            val2 = row2[idx2] if idx2 < len(row2) else ""

            if val1 != val2:
                changes["indices"].append(idx1)  # Use first file's index for reporting
                changes["column_names"].append(col_name)

                # Determine change type
                change_type = _classify_value_change(val1, val2)
                changes["change_types"].append(change_type)

                change_info = {
                    "column_index": idx1,
                    "column_name": col_name,
                    "change_type": change_type,
                    "value1": val1,
                    "value2": val2,
                    "value1_type": _get_value_type(val1),
                    "value2_type": _get_value_type(val2),
                    "value1_length": len(val1),
                    "value2_length": len(val2),
                }

                changes["field_changes"].append(change_info)
    else:
        # Fall back to index-based comparison when column structure differs
        max_len = max(len(row1), len(row2))

        for i in range(max_len):
            val1 = row1[i] if i < len(row1) else ""
            val2 = row2[i] if i < len(row2) else ""

            if val1 != val2:
                changes["indices"].append(i)

                # Determine change type
                change_type = _classify_value_change(val1, val2)
                changes["change_types"].append(change_type)

                # Get column name if available
                col_name = None
                if i < len(headers1):
                    col_name = headers1[i]
                elif i < len(headers2):
                    col_name = headers2[i]

                if col_name:
                    changes["column_names"].append(col_name)

                change_info = {
                    "column_index": i,
                    "column_name": col_name,
                    "change_type": change_type,
                    "value1": val1,
                    "value2": val2,
                    "value1_type": _get_value_type(val1),
                    "value2_type": _get_value_type(val2),
                    "value1_length": len(val1),
                    "value2_length": len(val2),
                }

                changes["field_changes"].append(change_info)

    return changes


def _extract_column_names_from_changes(changes: Dict[str, Any]) -> List[str]:
    """Extract column names from field changes, falling back to indices if names unavailable."""
    column_names = []
    for field_change in changes["field_changes"]:
        col_name = field_change.get("column_name")
        if col_name:
            column_names.append(col_name)
        else:
            # Fall back to index if column name is not available
            col_index = field_change.get("column_index", "unknown")
            column_names.append(f"column_{col_index}")
    return column_names


def _extract_header_names_from_changes(
    changes: Dict[str, Any], headers1: List[str], headers2: List[str]
) -> List[str]:
    """Extract header names from header changes, falling back to indices if names unavailable."""
    header_names = []
    for index in changes["indices"]:
        # Try to get header name from either row
        header_name = None
        if index < len(headers1):
            header_name = headers1[index]
        elif index < len(headers2):
            header_name = headers2[index]

        if header_name:
            header_names.append(header_name)
        else:
            header_names.append(f"column_{index}")
    return header_names


def _get_row_change_summary(changes: Dict[str, Any]) -> str:
    """Generate a human-readable summary of row changes."""
    if not changes["indices"]:
        return "no changes detected"

    change_counts: Dict[str, int] = {}
    for change_type in changes["change_types"]:
        change_counts[change_type] = change_counts.get(change_type, 0) + 1

    parts = []
    for change_type, count in change_counts.items():
        if count == 1:
            parts.append(f"1 {change_type} change")
        else:
            parts.append(f"{count} {change_type} changes")

    col_summary = f"in {len(changes['indices'])} column(s)"

    return f"{', '.join(parts)} {col_summary}"


def _analyze_column_differences(row1: List[str], row2: List[str]) -> tuple[int, int]:
    """Analyze column count differences between two rows."""
    len1, len2 = len(row1), len(row2)
    if len1 > len2:
        return len1 - len2, 0  # missing in run2, extra in run1
    elif len2 > len1:
        return 0, len2 - len1  # missing in run1, extra in run2
    else:
        return 0, 0


def _classify_value_change(val1: str, val2: str) -> str:
    """Classify the type of change between two values."""
    if not val1 and val2:
        return "added"
    elif val1 and not val2:
        return "removed"
    elif _is_numeric(val1) and _is_numeric(val2):
        return "numeric"
    elif val1.lower() in {
        "success",
        "failed",
        "complete",
        "incomplete",
    } or val2.lower() in {"success", "failed", "complete", "incomplete"}:
        return "outcome"
    else:
        return "text"


def _get_value_type(value: str) -> str:
    """Get the type classification of a value."""
    if not value:
        return "empty"
    elif _is_numeric(value):
        return "numeric"
    elif value.lower() in {
        "success",
        "failed",
        "complete",
        "incomplete",
        "true",
        "false",
    }:
        return "categorical"
    else:
        return "text"


def _is_numeric(value: str) -> bool:
    """Check if a string represents a numeric value."""
    try:
        float(value)
        return True
    except ValueError:
        return False


def _get_file_metadata(file_path: Path) -> Dict[str, Any]:
    """Get metadata information about a file."""
    try:
        stat = file_path.stat()
        return {
            "size": stat.st_size,
            "modified_time": datetime.fromtimestamp(stat.st_mtime).isoformat(),
            "exists": True,
        }
    except Exception:
        return {
            "size": 0,
            "modified_time": None,
            "exists": False,
        }


def _validate_column_names(csv_data: List[List[str]]) -> Optional[Dict[str, Any]]:
    """
    Validate that all column names in the CSV are unique.

    Args:
        csv_data: List of rows, where first row contains headers

    Returns:
        None if validation passes, otherwise a dict with error information
    """
    if not csv_data or not csv_data[0]:
        return None

    headers = csv_data[0]
    header_counts: Dict[str, int] = {}

    for i, header in enumerate(headers):
        if header in header_counts:
            return {
                "duplicate_header": header,
                "positions": [header_counts[header], i],
                "all_headers": headers,
            }
        header_counts[header] = i

    return None


def _create_column_mapping(headers: List[str]) -> Dict[str, int]:
    """
    Create a mapping from column names to their indices.

    Args:
        headers: List of column names

    Returns:
        Dictionary mapping column names to indices
    """
    return {header: i for i, header in enumerate(headers)}


def _get_column_value_by_name(
    row: List[str], column_name: str, column_map: Dict[str, int]
) -> str:
    """
    Get a value from a row by column name.

    Args:
        row: List of values in the row
        column_name: Name of the column to retrieve
        column_map: Mapping from column names to indices

    Returns:
        The value at the specified column, or empty string if not found
    """
    if column_name not in column_map:
        return ""

    col_idx = column_map[column_name]
    return row[col_idx] if col_idx < len(row) else ""


def _compare_column_orders(
    headers1: List[str], headers2: List[str]
) -> Optional[Dict[str, Any]]:
    """
    Compare the order of columns between two header lists.

    Args:
        headers1: First set of headers
        headers2: Second set of headers

    Returns:
        None if orders match, otherwise a dict with order difference information
    """
    # Only compare if both have the same columns
    set1 = set(headers1)
    set2 = set(headers2)

    if set1 != set2:
        return None  # Different columns, not just different order

    if headers1 == headers2:
        return None  # Same order

    # Find positions of each column in both files
    pos1 = {col: i for i, col in enumerate(headers1)}
    pos2 = {col: i for i, col in enumerate(headers2)}

    reordered_columns = []
    for col in headers1:
        if pos1[col] != pos2[col]:
            reordered_columns.append(
                {"column": col, "position_run1": pos1[col], "position_run2": pos2[col]}
            )

    return {
        "reordered_columns": reordered_columns,
        "headers_run1": headers1,
        "headers_run2": headers2,
    }


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
        file_size1 = file1.stat().st_size if file1.exists() else 0
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.FILE_READ_ERROR,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"First file is empty or failed to read (size: {file_size1} bytes)",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run1",
                    "file_path": str(file1),
                    "file_size": file_size1,
                },
            )
        )
        return discrepancies

    if not content2:
        file_size2 = file2.stat().st_size if file2.exists() else 0
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.FILE_READ_ERROR,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"Second file is empty or failed to read (size: {file_size2} bytes)",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run2",
                    "file_path": str(file2),
                    "file_size": file_size2,
                },
            )
        )
        return discrepancies

    # Get headers if available
    headers1 = content1[0] if content1 else []
    headers2 = content2[0] if content2 else []

    # Validate column names for duplicates
    dup_check1 = _validate_column_names(content1)
    if dup_check1:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.DUPLICATE_COLUMN_NAMES,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"Duplicate column name '{dup_check1['duplicate_header']}' at positions {dup_check1['positions']} in run1",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run1",
                    "duplicate_column": dup_check1["duplicate_header"],
                    "positions": dup_check1["positions"],
                },
                values={
                    "all_headers": dup_check1["all_headers"],
                    "duplicate_header": dup_check1["duplicate_header"],
                },
            )
        )

    dup_check2 = _validate_column_names(content2)
    if dup_check2:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.DUPLICATE_COLUMN_NAMES,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"Duplicate column name '{dup_check2['duplicate_header']}' at positions {dup_check2['positions']} in run2",
                location={
                    "file": filename,
                    "version": version,
                    "run": "run2",
                    "duplicate_column": dup_check2["duplicate_header"],
                    "positions": dup_check2["positions"],
                },
                values={
                    "all_headers": dup_check2["all_headers"],
                    "duplicate_header": dup_check2["duplicate_header"],
                },
            )
        )

    # Check for column order differences (only if no duplicate names and same columns)
    order_diff = _compare_column_orders(headers1, headers2)
    if order_diff:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.COLUMN_ORDER_DIFFERENCE,
                Severity.MEDIUM,
                Confidence.HIGH,
                f"Column order differs: {len(order_diff['reordered_columns'])} columns reordered",
                location={
                    "file": filename,
                    "version": version,
                    "run": "both",
                    "reordered_columns": [
                        col["column"] for col in order_diff["reordered_columns"]
                    ],
                },
                values={
                    "order_differences": order_diff,
                    "reordered_count": len(order_diff["reordered_columns"]),
                },
            )
        )

    # Create column mappings for name-based access (if no duplicates)
    column_map1 = _create_column_mapping(headers1) if not dup_check1 else {}
    column_map2 = _create_column_mapping(headers2) if not dup_check2 else {}

    # Discover index column for better comparison
    index_info = discover_index_column(content1, content2)

    # Log index column discovery results
    if index_info and index_info.get("index_column") is not None:
        logger.debug(
            f"Found index column for {filename}: {index_info['column_name']} (shared values: {index_info['shared_values']})"
        )
    else:
        logger.debug(
            f"No suitable index column found for {filename}: {index_info.get('reason', 'unknown') if index_info else 'no analysis'}"
        )

        # Add a low-severity discrepancy if no index column could be identified
        if (
            index_info
            and index_info.get("reason") in ["no_shared_values"]
            and index_info.get("candidate_columns", 0) > 0
        ):
            discrepancies.append(
                Discrepancy(
                    DiscrepancyType.NO_INDEX_COLUMN,
                    Severity.LOW,
                    Confidence.MEDIUM,
                    f"No suitable index column found: {index_info['candidate_columns']} unique columns but no shared values",
                    location={
                        "file": filename,
                        "version": version,
                        "run": "both",
                    },
                    values={
                        "index_analysis": index_info,
                        "candidate_columns": index_info.get("candidate_columns", 0),
                    },
                )
            )

    # Check if files have same number of rows
    if len(content1) != len(content2):
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.ROW_COUNT_DIFFERENCE,
                Severity.HIGH,
                Confidence.HIGH,
                f"Row count differs: {len(content1)} vs {len(content2)} (difference: {abs(len(content1) - len(content2))} rows)",
                location={"file": filename, "version": version, "run": "both"},
                values={
                    "run1_rows": len(content1),
                    "run2_rows": len(content2),
                    "row_difference": abs(len(content1) - len(content2)),
                },
            )
        )

    # Always check headers first, regardless of comparison method
    if headers1 != headers2:
        # Identify specific header differences
        changed_headers = _analyze_header_differences(headers1, headers2)
        header_details = _get_header_change_summary(changed_headers)
        changed_header_names = _extract_header_names_from_changes(
            changed_headers, headers1, headers2
        )

        discrepancies.append(
            Discrepancy(
                DiscrepancyType.HEADER_DIFFERENCE,
                Severity.CRITICAL,
                Confidence.HIGH,
                f"Header row differs: {header_details}",
                location={
                    "file": filename,
                    "version": version,
                    "run": "both",
                    "row": 1,
                    "changed_columns": changed_header_names,
                    "total_header_changes": len(changed_headers["indices"]),
                },
                values={
                    "header_changes": changed_headers,
                    "run1_header_count": len(headers1),
                    "run2_header_count": len(headers2),
                },
            )
        )

        # Check for column count differences in header
        if len(headers1) != len(headers2):
            missing_cols, extra_cols = _analyze_column_differences(headers1, headers2)
            discrepancies.append(
                Discrepancy(
                    DiscrepancyType.COLUMN_COUNT_DIFFERENCE,
                    Severity.CRITICAL,
                    Confidence.HIGH,
                    f"Header column count differs: {len(headers1)} vs {len(headers2)} ({missing_cols} missing, {extra_cols} extra in run2)",
                    location={
                        "file": filename,
                        "version": version,
                        "run": "both",
                        "row": 1,
                    },
                    values={
                        "run1_columns": len(headers1),
                        "run2_columns": len(headers2),
                        "columns_missing_in_run2": missing_cols,
                        "columns_extra_in_run2": extra_cols,
                    },
                )
            )

    # Use index-column-based comparison if available, otherwise fall back to position-based
    if (
        index_info
        and index_info.get("index_column") is not None
        and not dup_check1
        and not dup_check2
    ):
        # Use index-column-based row comparison (skip header row since already compared)
        index_column_name = index_info["column_name"]
        logger.debug(
            f"Using index-column-based comparison with column: {index_column_name}"
        )

        index_based_discrepancies = _compare_rows_by_index_column(
            content1,
            content2,
            index_column_name,
            headers1,
            headers2,
            column_map1,
            column_map2,
            version,
            filename,
        )
        discrepancies.extend(index_based_discrepancies)
    else:
        # Fall back to position-based row comparison (skip header row since already compared)
        logger.debug("Using position-based row comparison")

        # Compare data rows (skip header at index 0)
        max_rows = max(len(content1), len(content2))
        for i in range(1, max_rows):  # Start from 1 to skip header row
            row1 = content1[i] if i < len(content1) else []
            row2 = content2[i] if i < len(content2) else []

            if row1 != row2:
                # Analyze which specific columns differ in data rows
                column_differences = _analyze_row_differences(
                    row1, row2, headers1, headers2, column_map1, column_map2
                )

                # Determine severity based on content analysis
                severity = _determine_row_difference_severity(
                    row1, row2, i, column_differences
                )
                confidence = _determine_row_difference_confidence(row1, row2)

                change_summary = _get_row_change_summary(column_differences)
                changed_column_names = _extract_column_names_from_changes(
                    column_differences
                )

                discrepancies.append(
                    Discrepancy(
                        DiscrepancyType.ROW_DIFFERENCE,
                        severity,
                        confidence,
                        f"Row {i + 1} differs: {change_summary}",
                        location={
                            "file": filename,
                            "version": version,
                            "run": "both",
                            "row": i + 1,
                            "changed_columns": changed_column_names,
                            "total_field_changes": len(column_differences["indices"]),
                        },
                        values={
                            "field_changes": column_differences,
                            "change_types": column_differences["change_types"],
                        },
                    )
                )

                # Check for column count differences in data rows
                if len(row1) != len(row2):
                    missing_cols, extra_cols = _analyze_column_differences(row1, row2)
                    discrepancies.append(
                        Discrepancy(
                            DiscrepancyType.COLUMN_COUNT_DIFFERENCE,
                            Severity.HIGH,
                            Confidence.HIGH,
                            f"Row {i + 1} column count differs: {len(row1)} vs {len(row2)} ({missing_cols} missing, {extra_cols} extra in run2)",
                            location={
                                "file": filename,
                                "version": version,
                                "run": "both",
                                "row": i + 1,
                                "column_difference": abs(len(row1) - len(row2)),
                            },
                            values={
                                "run1_columns": len(row1),
                                "run2_columns": len(row2),
                                "columns_missing_in_run2": missing_cols,
                                "columns_extra_in_run2": extra_cols,
                            },
                        )
                    )

    return discrepancies


def _determine_row_difference_severity(
    row1: List[str], row2: List[str], row_index: int, column_differences: Dict[str, Any]
) -> Severity:
    """Determine severity of row difference based on content analysis."""
    # Check for outcome changes (highest priority)
    for change in column_differences.get("field_changes", []):
        if change["change_type"] == "outcome":
            return Severity.HIGH

    # Check for numeric changes
    numeric_changes = sum(
        1
        for change in column_differences.get("field_changes", [])
        if change["change_type"] == "numeric"
    )
    if numeric_changes > 0:
        return Severity.MEDIUM

    # Check if many fields changed
    total_changes = len(column_differences.get("indices", []))
    if total_changes > 3:
        return Severity.MEDIUM

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


def _build_index_row_mapping(
    csv_data: List[List[str]], index_column_name: str
) -> Dict[str, Dict[str, Any]]:
    """
    Build a mapping from index column values to row data and positions.

    Args:
        csv_data: List of rows, where first row contains headers
        index_column_name: Name of the column to use as index

    Returns:
        Dictionary mapping index values to {"row": row_data, "position": row_position}
    """
    if not csv_data or len(csv_data) < 2:  # Need at least header + 1 data row
        return {}

    headers = csv_data[0]
    if index_column_name not in headers:
        return {}

    index_col_idx = headers.index(index_column_name)
    mapping = {}

    # Start from row 1 to skip header
    for row_idx in range(1, len(csv_data)):
        row = csv_data[row_idx]
        if index_col_idx < len(row):
            index_value = row[index_col_idx].strip()
            if index_value:  # Only map non-empty index values
                mapping[index_value] = {"row": row, "position": row_idx}

    return mapping


def _compare_rows_by_index_column(
    content1: List[List[str]],
    content2: List[List[str]],
    index_column_name: str,
    headers1: List[str],
    headers2: List[str],
    column_map1: Dict[str, int],
    column_map2: Dict[str, int],
    version: str,
    filename: str,
) -> List[Discrepancy]:
    """
    Compare CSV rows using index column matching instead of position-based comparison.

    Args:
        content1: First CSV data
        content2: Second CSV data
        index_column_name: Name of column to use for row matching
        headers1: Headers from first file
        headers2: Headers from second file
        column_map1: Column name to index mapping for first file
        column_map2: Column name to index mapping for second file
        version: Version being compared
        filename: Name of the file being compared

    Returns:
        List of discrepancies found during comparison
    """
    discrepancies = []

    # Build index-to-row mappings
    index_map1 = _build_index_row_mapping(content1, index_column_name)
    index_map2 = _build_index_row_mapping(content2, index_column_name)

    # Get all unique index values from both files
    all_index_values = set(index_map1.keys()) | set(index_map2.keys())

    # Track row position differences for order comparison
    position_differences = []

    for index_value in sorted(all_index_values):
        row1_info = index_map1.get(index_value)
        row2_info = index_map2.get(index_value)

        if row1_info and row2_info:
            # Row exists in both files - compare content and positions
            row1 = row1_info["row"]
            row2 = row2_info["row"]
            pos1 = row1_info["position"]
            pos2 = row2_info["position"]

            # Check if row positions differ (for order difference detection)
            if pos1 != pos2:
                position_differences.append(
                    {
                        "index_value": index_value,
                        "position_run1": pos1,
                        "position_run2": pos2,
                    }
                )

            # Compare row content if rows differ
            if row1 != row2:
                # Analyze which specific columns differ
                column_differences = _analyze_row_differences(
                    row1, row2, headers1, headers2, column_map1, column_map2
                )

                # Determine severity based on content analysis
                severity = _determine_row_difference_severity(
                    row1, row2, pos1, column_differences
                )
                confidence = _determine_row_difference_confidence(row1, row2)

                change_summary = _get_row_change_summary(column_differences)
                changed_column_names = _extract_column_names_from_changes(
                    column_differences
                )

                discrepancies.append(
                    Discrepancy(
                        DiscrepancyType.ROW_DIFFERENCE,
                        severity,
                        confidence,
                        f"Row with {index_column_name}='{index_value}' differs: {change_summary}",
                        location={
                            "file": filename,
                            "version": version,
                            "run": "both",
                            "row": pos1
                            + 1,  # Add backward compatibility field (1-based row number)
                            "index_column": index_column_name,
                            "index_value": index_value,
                            "position_run1": pos1,
                            "position_run2": pos2,
                            "changed_columns": changed_column_names,
                            "total_field_changes": len(column_differences["indices"]),
                        },
                        values={
                            "field_changes": column_differences,
                            "change_types": column_differences["change_types"],
                        },
                    )
                )

                # Check for column count differences in data rows
                if len(row1) != len(row2):
                    missing_cols, extra_cols = _analyze_column_differences(row1, row2)
                    discrepancies.append(
                        Discrepancy(
                            DiscrepancyType.COLUMN_COUNT_DIFFERENCE,
                            Severity.HIGH,
                            Confidence.HIGH,
                            f"Row with {index_column_name}='{index_value}' column count differs: {len(row1)} vs {len(row2)} ({missing_cols} missing, {extra_cols} extra in run2)",
                            location={
                                "file": filename,
                                "version": version,
                                "run": "both",
                                "index_column": index_column_name,
                                "index_value": index_value,
                                "position_run1": pos1,
                                "position_run2": pos2,
                                "column_difference": abs(len(row1) - len(row2)),
                            },
                            values={
                                "run1_columns": len(row1),
                                "run2_columns": len(row2),
                                "columns_missing_in_run2": missing_cols,
                                "columns_extra_in_run2": extra_cols,
                            },
                        )
                    )

        elif row1_info and not row2_info:
            # Row exists only in first file
            row1 = row1_info["row"]
            pos1 = row1_info["position"]

            discrepancies.append(
                Discrepancy(
                    DiscrepancyType.MISSING_ROW,
                    Severity.HIGH,
                    Confidence.HIGH,
                    f"Row with {index_column_name}='{index_value}' missing in run2",
                    location={
                        "file": filename,
                        "version": version,
                        "run": "run2",
                        "index_column": index_column_name,
                        "index_value": index_value,
                        "position_run1": pos1,
                        "missing_from": "run2",
                        "present_in": "run1",
                    },
                    values={
                        "missing_row_data": row1,
                    },
                )
            )

        elif row2_info and not row1_info:
            # Row exists only in second file
            row2 = row2_info["row"]
            pos2 = row2_info["position"]

            discrepancies.append(
                Discrepancy(
                    DiscrepancyType.EXTRA_ROW,
                    Severity.HIGH,
                    Confidence.HIGH,
                    f"Row with {index_column_name}='{index_value}' extra in run2",
                    location={
                        "file": filename,
                        "version": version,
                        "run": "run2",
                        "index_column": index_column_name,
                        "index_value": index_value,
                        "position_run2": pos2,
                        "missing_from": "run1",
                        "present_in": "run2",
                    },
                    values={
                        "extra_row_data": row2,
                    },
                )
            )

    # Report row order differences if rows exist in both files but in different positions
    if position_differences:
        discrepancies.append(
            Discrepancy(
                DiscrepancyType.ROW_ORDER_DIFFERENCE,
                Severity.LOW,
                Confidence.HIGH,
                f"Row order differs: {len(position_differences)} rows in different positions",
                location={
                    "file": filename,
                    "version": version,
                    "run": "both",
                    "index_column": index_column_name,
                    "reordered_rows": [
                        diff["index_value"] for diff in position_differences
                    ],
                },
                values={
                    "position_differences": position_differences,
                    "reordered_count": len(position_differences),
                },
            )
        )

    return discrepancies


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
            # Get file info from run2
            file2_info = _get_file_metadata(csv_files2[csv_name])
            report.add_discrepancy(
                version,
                csv_name,
                Discrepancy(
                    DiscrepancyType.MISSING_FILE,
                    Severity.CRITICAL,
                    Confidence.HIGH,
                    f"File missing in run1: {csv_name} (run2 size: {file2_info['size']} bytes)",
                    location={
                        "file": csv_name,
                        "version": version,
                        "run": "run1",
                        "missing_from": "run1",
                        "present_in": "run2",
                    },
                    values={
                        "file_info_run2": file2_info,
                    },
                ),
            )
            continue

        if csv_name not in csv_files2:
            # Get file info from run1
            file1_info = _get_file_metadata(csv_files1[csv_name])
            report.add_discrepancy(
                version,
                csv_name,
                Discrepancy(
                    DiscrepancyType.MISSING_FILE,
                    Severity.CRITICAL,
                    Confidence.HIGH,
                    f"File missing in run2: {csv_name} (run1 size: {file1_info['size']} bytes)",
                    location={
                        "file": csv_name,
                        "version": version,
                        "run": "run2",
                        "missing_from": "run2",
                        "present_in": "run1",
                    },
                    values={
                        "file_info_run1": file1_info,
                    },
                ),
            )
            continue

        # Compare the files
        discrepancies = compare_csv_contents(
            csv_files1[csv_name], csv_files2[csv_name], version, csv_name
        )

        if not discrepancies:
            # For identical files, we still want to discover and store index column info
            content1 = read_csv_file(csv_files1[csv_name])
            content2 = read_csv_file(csv_files2[csv_name])
            index_info = (
                discover_index_column(content1, content2)
                if content1 and content2
                else None
            )
            report.mark_file_identical(version, csv_name, index_info)
        else:
            # Set index column info for files with differences
            content1 = read_csv_file(csv_files1[csv_name])
            content2 = read_csv_file(csv_files2[csv_name])
            index_info = (
                discover_index_column(content1, content2)
                if content1 and content2
                else None
            )

            for discrepancy in discrepancies:
                report.add_discrepancy(version, csv_name, discrepancy)

            # Set the index column info after adding discrepancies
            report.set_file_index_info(version, csv_name, index_info)


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


def _find_unique_value_columns(csv_data: List[List[str]]) -> List[str]:
    """
    Find all columns that have unique values (no duplicates).

    Args:
        csv_data: List of rows, where each row is a list of values

    Returns:
        List of column names that contain only unique values
    """
    if not csv_data or len(csv_data) < 2:  # Need at least header + 1 data row
        return []

    unique_columns = []
    headers = csv_data[0] if csv_data else []
    num_columns = len(headers)

    for col_idx in range(num_columns):
        # Extract all values from this column (skip header row)
        column_values = []
        for row_idx in range(1, len(csv_data)):  # Skip header
            if col_idx < len(csv_data[row_idx]):
                value = csv_data[row_idx][col_idx].strip()
                if value:  # Only consider non-empty values
                    column_values.append(value)

        # Check if all values are unique
        if len(column_values) > 0 and len(column_values) == len(set(column_values)):
            col_name = (
                headers[col_idx] if col_idx < len(headers) else f"column_{col_idx}"
            )
            unique_columns.append(col_name)

    return unique_columns


def _count_shared_values(
    csv_data1: List[List[str]], csv_data2: List[List[str]], column_name: str
) -> int:
    """
    Count how many values are shared between two CSV files in a specific column.

    Args:
        csv_data1: First CSV data
        csv_data2: Second CSV data
        column_name: Column name to compare

    Returns:
        Number of shared non-empty values in the specified column
    """
    if not csv_data1 or not csv_data2:
        return 0

    headers1 = csv_data1[0] if csv_data1 else []
    headers2 = csv_data2[0] if csv_data2 else []

    # Find column indices by name
    col_idx1 = headers1.index(column_name) if column_name in headers1 else -1
    col_idx2 = headers2.index(column_name) if column_name in headers2 else -1

    if col_idx1 == -1 or col_idx2 == -1:
        return 0

    # Extract values from column in first file (skip header)
    values1 = set()
    for row_idx in range(1, len(csv_data1)):
        if col_idx1 < len(csv_data1[row_idx]):
            value = csv_data1[row_idx][col_idx1].strip()
            if value:
                values1.add(value)

    # Extract values from column in second file (skip header)
    values2 = set()
    for row_idx in range(1, len(csv_data2)):
        if col_idx2 < len(csv_data2[row_idx]):
            value = csv_data2[row_idx][col_idx2].strip()
            if value:
                values2.add(value)

    # Count intersection
    return len(values1 & values2)


def discover_index_column(
    csv_data1: List[List[str]], csv_data2: List[List[str]]
) -> Optional[Dict[str, Any]]:
    """
    Discover the best index column for comparing two CSV files.

    An index column must:
    1. Have all unique values within each file
    2. Have the highest number of values shared between the two files

    Args:
        csv_data1: First CSV data as list of rows
        csv_data2: Second CSV data as list of rows

    Returns:
        Dictionary with index column information, or None if no suitable column found
    """
    if not csv_data1 or not csv_data2:
        return None

    # Get headers from both files
    headers1 = csv_data1[0] if csv_data1 else []
    headers2 = csv_data2[0] if csv_data2 else []

    # Find columns with unique values in both files
    unique_cols1 = set(_find_unique_value_columns(csv_data1))
    unique_cols2 = set(_find_unique_value_columns(csv_data2))

    # Only consider columns that are unique in both files
    candidate_columns = unique_cols1 & unique_cols2

    if not candidate_columns:
        return None

    # Find the column with the most shared values
    best_column = None
    max_shared = 0
    column_analysis = {}

    for col_name in sorted(candidate_columns):
        shared_count = _count_shared_values(csv_data1, csv_data2, col_name)

        # Get column indices for analysis
        col_idx1 = headers1.index(col_name) if col_name in headers1 else -1
        col_idx2 = headers2.index(col_name) if col_name in headers2 else -1

        column_analysis[col_name] = {
            "column_name": col_name,
            "column_index_run1": col_idx1,
            "column_index_run2": col_idx2,
            "shared_values": shared_count,
            "unique_in_both": True,
        }

        if shared_count > max_shared:
            max_shared = shared_count
            best_column = col_name

    if best_column is None:
        return {
            "index_column": None,
            "column_name": None,
            "shared_values": 0,
            "candidate_columns": len(candidate_columns),
            "column_analysis": column_analysis,
            "reason": "no_shared_values",
        }

    # Get the best column info
    best_col_info = column_analysis[best_column]

    return {
        "index_column": best_col_info[
            "column_index_run1"
        ],  # Return index for backward compatibility
        "column_name": best_column,
        "shared_values": max_shared,
        "candidate_columns": len(candidate_columns),
        "column_analysis": column_analysis,
        "reason": "success",
    }


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
