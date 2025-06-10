"""
Discrepancy classes for the compare_runs script.

This module defines the core classes used to represent and manage
discrepancies found during proviral analysis comparison. It uses
a class-based approach with subclassing instead of Enums for better
type safety and extensibility.

The module provides:
- Severity levels (Critical, High, Medium, Low) via the Severity class
- Confidence levels (High, Medium, Low) via the Confidence class
- Main Discrepancy and ComparisonReport classes for managing findings
"""

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from enum import Enum


# Severity levels for discrepancies
class Severity(Enum):
    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


# Confidence levels for discrepancies
class Confidence(Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


def _trim_value_for_display(value: str, max_length: int = 20) -> str:
    """
    Trim a value for display purposes if it exceeds the maximum length.

    For values longer than max_length, truncates with "..." in the middle
    to preserve both the beginning and end of the value.

    Args:
        value: The string value to potentially trim
        max_length: Maximum allowed length (default: 20)

    Returns:
        The original value if <= max_length, otherwise a trimmed version

    Examples:
        >>> _trim_value_for_display("short")
        'short'
        >>> _trim_value_for_display("this_is_a_very_long_value_that_exceeds_twenty_characters")
        'this_is_...characters'
    """
    if not isinstance(value, str):
        value = str(value)

    if len(value) <= max_length:
        return value

    # Calculate prefix and suffix lengths
    # We need room for "..." (3 chars) in the middle
    available_chars = max_length - 3
    prefix_length = available_chars // 2
    suffix_length = available_chars - prefix_length

    # Ensure we don't get negative lengths for very small max_length
    if prefix_length < 0 or suffix_length < 0:
        return value[:max_length]

    prefix = value[:prefix_length]
    suffix = value[-suffix_length:] if suffix_length > 0 else ""

    return f"{prefix}...{suffix}"


def _trim_data_recursively(data: Any, max_length: int = 20) -> Any:
    """
    Recursively trim string values in a data structure (dict, list, or string).
    """
    if isinstance(data, str):
        return _trim_value_for_display(data, max_length)
    elif isinstance(data, list):
        return [_trim_data_recursively(item, max_length) for item in data]
    elif isinstance(data, dict):
        return {
            key: _trim_data_recursively(value, max_length)
            for key, value in data.items()
        }
    return data


# Base class for all discrepancy types
class DiscrepancyBase:
    """Base class for all specific discrepancy types."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
    ):
        self.severity = severity
        self.confidence = confidence
        self.description = description
        # Common location fields
        self.file = file
        self.version = version
        self.run = run

    def to_dict(self) -> Dict[str, Any]:
        """Convert discrepancy to dictionary for JSON output."""
        location_dict: Dict[str, Any] = {
            "file": self.file,
            "version": self.version,
            "run": self.run,
        }
        values_dict: Dict[str, Any] = {}

        # Add type-specific location and values fields
        self._add_location_fields(location_dict)
        self._add_values_fields(values_dict)
        values = _trim_data_recursively(values_dict)

        return {
            "type": self.__class__.__name__,
            "severity": self.severity.value,
            "confidence": self.confidence.value,
            "description": self.description,
            "location": _trim_data_recursively(location_dict),
            **values,
        }

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        """Add type-specific location fields. Override in subclasses."""
        pass

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        """Add type-specific values fields. Override in subclasses."""
        pass


class MissingFileDiscrepancy(DiscrepancyBase):
    """Represents a missing file discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        missing_from: str,
        present_in: str,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.missing_from = missing_from
        self.present_in = present_in

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["missing_from"] = self.missing_from

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["missing_from"] = self.missing_from
        values_dict["present_in"] = self.present_in


class MissingDirectoryDiscrepancy(DiscrepancyBase):
    """Represents a missing directory discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        missing_from: str,
        present_in: str,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.missing_from = missing_from
        self.present_in = present_in

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["missing_from"] = self.missing_from

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["missing_from"] = self.missing_from
        values_dict["present_in"] = self.present_in


class HeaderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a header difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        changed_headers: List[str],
        total_header_changes: int,
        header_changes: Dict[str, Any],
        run1_header_count: int,
        run2_header_count: int,
        row: int = 1,  # Headers are typically row 1
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.changed_headers = changed_headers
        self.total_header_changes = total_header_changes
        self.header_changes = header_changes
        self.run1_header_count = run1_header_count
        self.run2_header_count = run2_header_count
        self.row = row

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["row"] = self.row
        location_dict["changed_headers"] = self.changed_headers
        location_dict["total_header_changes"] = self.total_header_changes

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["header_changes"] = self.header_changes
        values_dict["run1_header_count"] = self.run1_header_count
        values_dict["run2_header_count"] = self.run2_header_count


class RowDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        row: int,
        index_column: str,
        index_value: str,
        position_run1: int,
        position_run2: int,
        changed_columns: List[str],
        total_field_changes: int,
        field_changes: Dict[str, Any],
        change_types: List[str],
    ):
        super().__init__(
            severity, confidence, description, file, version, run
        )
        self.row = row
        self.index_column = index_column
        self.index_value = index_value
        self.position_run1 = position_run1
        self.position_run2 = position_run2
        self.changed_columns = changed_columns
        self.total_field_changes = total_field_changes
        self.field_changes = field_changes
        self.change_types = change_types

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["row"] = self.row
        location_dict["index_column"] = self.index_column
        location_dict["index_value"] = self.index_value
        location_dict["position_run1"] = self.position_run1
        location_dict["position_run2"] = self.position_run2
        location_dict["changed_columns"] = self.changed_columns
        location_dict["total_field_changes"] = self.total_field_changes

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["field_changes"] = self.field_changes
        values_dict["change_types"] = self.change_types


class RowCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row count difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        run1_rows: int,
        run2_rows: int,
        row_difference: int,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.run1_rows = run1_rows
        self.run2_rows = run2_rows
        self.row_difference = row_difference

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["row_difference"] = self.row_difference

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["run1_rows"] = self.run1_rows
        values_dict["run2_rows"] = self.run2_rows
        values_dict["row_difference"] = self.row_difference


class ColumnCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column count difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        run1_columns: int,
        run2_columns: int,
        columns_missing_in_run2: int,
        columns_extra_in_run2: int,
        row: Optional[int] = None,
        index_column: Optional[str] = None,
        index_value: Optional[str] = None,
        position_run1: Optional[int] = None,
        position_run2: Optional[int] = None,
        column_difference: Optional[int] = None,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.run1_columns = run1_columns
        self.run2_columns = run2_columns
        self.columns_missing_in_run2 = columns_missing_in_run2
        self.columns_extra_in_run2 = columns_extra_in_run2
        # Optional fields for data row column count differences
        self.row = row
        self.index_column = index_column
        self.index_value = index_value
        self.position_run1 = position_run1
        self.position_run2 = position_run2
        self.column_difference = column_difference

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        if self.row is not None:
            location_dict["row"] = self.row
        if self.index_column is not None:
            location_dict["index_column"] = self.index_column
        if self.index_value is not None:
            location_dict["index_value"] = self.index_value
        if self.position_run1 is not None:
            location_dict["position_run1"] = self.position_run1
        if self.position_run2 is not None:
            location_dict["position_run2"] = self.position_run2
        if self.column_difference is not None:
            location_dict["column_difference"] = self.column_difference

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["run1_columns"] = self.run1_columns
        values_dict["run2_columns"] = self.run2_columns
        values_dict["columns_missing_in_run2"] = self.columns_missing_in_run2
        values_dict["columns_extra_in_run2"] = self.columns_extra_in_run2


class FileReadErrorDiscrepancy(DiscrepancyBase):
    """Represents a file read error discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        file_path: str,
        file_size: Optional[int] = None,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.file_path = file_path
        self.file_size = file_size

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["file_path"] = self.file_path
        if self.file_size is not None:
            location_dict["file_size"] = self.file_size

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["file_path"] = self.file_path
        if self.file_size is not None:
            values_dict["file_size"] = self.file_size


class EmptyFileDiscrepancy(DiscrepancyBase):
    """Represents an empty file discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        file_path: str,
        file_size: int,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.file_path = file_path
        self.file_size = file_size

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["file_path"] = self.file_path
        location_dict["file_size"] = self.file_size

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["file_path"] = self.file_path
        values_dict["file_size"] = self.file_size


class NoIndexColumnDiscrepancy(DiscrepancyBase):
    """Represents a no index column discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        reason: str,
        index_analysis: Optional[Dict[str, Any]],
        has_duplicate_columns_run1: bool,
        has_duplicate_columns_run2: bool,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.reason = reason
        self.index_analysis = index_analysis
        self.has_duplicate_columns_run1 = has_duplicate_columns_run1
        self.has_duplicate_columns_run2 = has_duplicate_columns_run2

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["reason"] = self.reason

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["index_analysis"] = self.index_analysis
        values_dict["has_duplicate_columns_run1"] = self.has_duplicate_columns_run1
        values_dict["has_duplicate_columns_run2"] = self.has_duplicate_columns_run2
        values_dict["reason"] = self.reason


class DuplicateColumnNamesDiscrepancy(DiscrepancyBase):
    """Represents a duplicate column names discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        duplicate_column: str,
        positions: List[int],
        all_headers: List[str],
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.duplicate_column = duplicate_column
        self.positions = positions
        self.all_headers = all_headers

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["duplicate_column"] = self.duplicate_column
        location_dict["positions"] = self.positions

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["all_headers"] = self.all_headers
        values_dict["duplicate_header"] = self.duplicate_column


class ColumnOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column order difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        reordered_columns: List[str],
        order_differences: Dict[str, Any],
        reordered_count: int,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.reordered_columns = reordered_columns
        self.order_differences = order_differences
        self.reordered_count = reordered_count

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["reordered_columns"] = self.reordered_columns

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["order_differences"] = self.order_differences
        values_dict["reordered_count"] = self.reordered_count


class RowOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row order difference discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        index_column: str,
        reordered_rows: List[str],
        position_differences: List[Dict[str, Any]],
        reordered_count: int,
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.index_column = index_column
        self.reordered_rows = reordered_rows
        self.position_differences = position_differences
        self.reordered_count = reordered_count

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["index_column"] = self.index_column
        location_dict["reordered_rows"] = self.reordered_rows

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["position_differences"] = self.position_differences
        values_dict["reordered_count"] = self.reordered_count


class MissingRowDiscrepancy(DiscrepancyBase):
    """Represents a missing row discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        index_column: str,
        index_value: str,
        position_run1: int,
        missing_from: str,
        present_in: str,
        missing_row_data: List[str],
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.index_column = index_column
        self.index_value = index_value
        self.position_run1 = position_run1
        self.missing_from = missing_from
        self.present_in = present_in
        self.missing_row_data = missing_row_data

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["index_column"] = self.index_column
        location_dict["index_value"] = self.index_value
        location_dict["position_run1"] = self.position_run1
        location_dict["missing_from"] = self.missing_from
        location_dict["present_in"] = self.present_in

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["missing_row_data"] = self.missing_row_data


class ExtraRowDiscrepancy(DiscrepancyBase):
    """Represents an extra row discrepancy."""

    def __init__(
        self,
        severity: Severity,
        confidence: Confidence,
        description: str,
        file: str,
        version: str,
        run: str,
        index_column: str,
        index_value: str,
        position_run2: int,
        missing_from: str,
        present_in: str,
        extra_row_data: List[str],
    ):
        super().__init__(severity, confidence, description, file, version, run)
        self.index_column = index_column
        self.index_value = index_value
        self.position_run2 = position_run2
        self.missing_from = missing_from
        self.present_in = present_in
        self.extra_row_data = extra_row_data

    def _add_location_fields(self, location_dict: Dict[str, Any]) -> None:
        location_dict["index_column"] = self.index_column
        location_dict["index_value"] = self.index_value
        location_dict["position_run2"] = self.position_run2
        location_dict["missing_from"] = self.missing_from
        location_dict["present_in"] = self.present_in

    def _add_values_fields(self, values_dict: Dict[str, Any]) -> None:
        values_dict["extra_row_data"] = self.extra_row_data


# Union type for all specific discrepancy classes
Discrepancy = Union[
    MissingFileDiscrepancy,
    MissingDirectoryDiscrepancy,
    HeaderDifferenceDiscrepancy,
    RowDifferenceDiscrepancy,
    RowCountDifferenceDiscrepancy,
    ColumnCountDifferenceDiscrepancy,
    FileReadErrorDiscrepancy,
    EmptyFileDiscrepancy,
    NoIndexColumnDiscrepancy,
    DuplicateColumnNamesDiscrepancy,
    ColumnOrderDifferenceDiscrepancy,
    RowOrderDifferenceDiscrepancy,
    MissingRowDiscrepancy,
    ExtraRowDiscrepancy,
]


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
