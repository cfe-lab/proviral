"""
Discrepancy classes for the compare_runs script.

This module defines the core classes used to represent and manage
discrepancies found during proviral analysis comparison. It uses
a dataclass-based approach with automatic serialization for better
type safety and cleaner code.

The module provides:
- Severity levels (Critical, High, Medium, Low) via the Severity class
- Confidence levels (High, Medium, Low) via the Confidence class
- Location field decorator for automatic categorization
- Main Discrepancy and ComparisonReport classes for managing findings
"""

import json
from collections import defaultdict
from dataclasses import dataclass, field, fields
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from enum import Enum


# Location field marker for dataclass fields that should go in location subobject
def location_field(**kwargs):
    """Mark a dataclass field as a location field."""
    kwargs["metadata"] = {"is_location": True}
    return field(**kwargs)


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
@dataclass(frozen=True)
class DiscrepancyBase:
    """Base class for all specific discrepancy types."""

    # Standard fields
    severity: Severity
    confidence: Confidence
    description: str
    file: str
    version: str
    run: str

    def to_dict(self) -> Dict[str, Any]:
        """Convert discrepancy to dictionary for JSON output."""
        location_dict: Dict[str, Any] = {}
        top_level_dict: Dict[str, Any] = {}

        # Standard location fields from base class
        base_location_fields = {"file", "version", "run"}

        # Additional location fields for specific discrepancy types
        column_count_location_fields = {
            "row",
            "index_column",
            "index_value",
            "position_run1",
            "position_run2",
            "column_difference",
        }

        # Traverse all dataclass fields
        for field_info in fields(self):
            value = getattr(self, field_info.name)

            # Skip None values for optional fields
            if value is None:
                continue

            # Check if this is a location field
            is_location = (
                field_info.metadata.get("is_location", False)
                or field_info.name in base_location_fields
                or (
                    isinstance(self, ColumnCountDifferenceDiscrepancy)
                    and field_info.name in column_count_location_fields
                )
            )

            if field_info.name in ["severity", "confidence"]:
                # Handle enum fields
                top_level_dict[field_info.name] = value.value
            elif field_info.name == "description":
                # Don't trim description
                top_level_dict[field_info.name] = value
            elif is_location:
                # Add to location dict and trim
                location_dict[field_info.name] = value
            else:
                # Add to top level and trim
                top_level_dict[field_info.name] = value

        # Special handling for MissingFileDiscrepancy and MissingDirectoryDiscrepancy
        # missing_from should appear in both location and top level
        if isinstance(self, (MissingFileDiscrepancy, MissingDirectoryDiscrepancy)):
            if hasattr(self, "missing_from"):
                location_dict["missing_from"] = self.missing_from

        # Build final result
        result_dict = {
            "type": self.__class__.__name__,
            "location": _trim_data_recursively(location_dict),
            **{k: v for k, v in top_level_dict.items() if k != "description"},
        }

        # Add description untrimmed
        if "description" in top_level_dict:
            result_dict["description"] = top_level_dict["description"]

        # Trim everything except description and location
        return _trim_data_recursively(
            {
                k: v
                for k, v in result_dict.items()
                if k not in ["description", "location"]
            }
        ) | {
            "description": result_dict.get("description", ""),
            "location": result_dict["location"],
        }


@dataclass(frozen=True)
class MissingFileDiscrepancy(DiscrepancyBase):
    """Represents a missing file discrepancy."""

    missing_from: str
    present_in: str


@dataclass(frozen=True)
class MissingDirectoryDiscrepancy(DiscrepancyBase):
    """Represents a missing directory discrepancy."""

    missing_from: str
    present_in: str


@dataclass(frozen=True)
class HeaderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a header difference discrepancy."""

    header_changes: Dict[str, Any]
    run1_header_count: int
    run2_header_count: int
    row: int = location_field()
    changed_headers: List[str] = location_field()
    total_header_changes: int = location_field()


@dataclass(frozen=True)
class RowDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row difference discrepancy."""

    field_changes: Dict[str, Any]
    change_types: List[str]
    row: int = location_field()
    index_column: str = location_field()
    index_value: str = location_field()
    position_run1: int = location_field()
    position_run2: int = location_field()
    changed_columns: List[str] = location_field()
    total_field_changes: int = location_field()


@dataclass(frozen=True)
class RowCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row count difference discrepancy."""

    run1_rows: int
    run2_rows: int
    row_difference: int = location_field()


@dataclass(frozen=True)
class ColumnCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column count difference discrepancy."""

    run1_columns: int
    run2_columns: int
    columns_missing_in_run2: int
    columns_extra_in_run2: int
    # Optional fields for data row column count differences
    row: Optional[int]
    index_column: Optional[str]
    index_value: Optional[str]
    position_run1: Optional[int]
    position_run2: Optional[int]
    column_difference: Optional[int]


@dataclass(frozen=True)
class FileReadErrorDiscrepancy(DiscrepancyBase):
    """Represents a file read error discrepancy."""

    file_path: str
    file_size: Optional[int]

    def to_dict(self) -> Dict[str, Any]:
        """Override to add file_path to location."""
        result = super().to_dict()
        # Add file_path to location
        result["location"]["file_path"] = self.file_path
        if self.file_size is not None:
            result["location"]["file_size"] = self.file_size
        return result


@dataclass(frozen=True)
class EmptyFileDiscrepancy(DiscrepancyBase):
    """Represents an empty file discrepancy."""

    file_path: str
    file_size: int

    def to_dict(self) -> Dict[str, Any]:
        """Override to add file_path and file_size to location."""
        result = super().to_dict()
        # Add to location
        result["location"]["file_path"] = self.file_path
        result["location"]["file_size"] = self.file_size
        return result


@dataclass(frozen=True)
class NoIndexColumnDiscrepancy(DiscrepancyBase):
    """Represents a no index column discrepancy."""

    reason: str
    index_analysis: Optional[Dict[str, Any]]
    has_duplicate_columns_run1: bool
    has_duplicate_columns_run2: bool

    def to_dict(self) -> Dict[str, Any]:
        """Override to add reason to location."""
        result = super().to_dict()
        # Add reason to location
        result["location"]["reason"] = self.reason
        return result


@dataclass(frozen=True)
class DuplicateColumnNamesDiscrepancy(DiscrepancyBase):
    """Represents a duplicate column names discrepancy."""

    all_headers: List[str]
    duplicate_header: str  # alias for duplicate_column
    duplicate_column: str  # Original field name
    positions: List[int]

    def to_dict(self) -> Dict[str, Any]:
        """Override to add duplicate_column to location and handle aliasing."""
        result = super().to_dict()
        # Add to location
        result["location"]["duplicate_column"] = self.duplicate_column
        result["location"]["positions"] = self.positions
        # Set alias field for backward compatibility
        result["duplicate_header"] = self.duplicate_column
        return result


@dataclass(frozen=True)
class ColumnOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column order difference discrepancy."""

    order_differences: Dict[str, Any]
    reordered_count: int
    reordered_columns: List[str] = location_field()


@dataclass(frozen=True)
class RowOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row order difference discrepancy."""

    position_differences: List[Dict[str, Any]]
    reordered_count: int
    index_column: str = location_field()
    reordered_rows: List[str] = location_field()


@dataclass(frozen=True)
class MissingRowDiscrepancy(DiscrepancyBase):
    """Represents a missing row discrepancy."""

    missing_row_data: List[str]
    index_column: str = location_field()
    index_value: str = location_field()
    position_run1: int = location_field()
    missing_from: str = location_field()
    present_in: str = location_field()


@dataclass(frozen=True)
class ExtraRowDiscrepancy(DiscrepancyBase):
    """Represents an extra row discrepancy."""

    extra_row_data: List[str]
    index_column: str = location_field()
    index_value: str = location_field()
    position_run2: int = location_field()
    missing_from: str = location_field()
    present_in: str = location_field()


@dataclass(frozen=True)
class FieldChangeDiscrepancy(DiscrepancyBase):
    """Represents a single field change within a row."""

    change_type: str
    value1: str
    value2: str
    value1_type: str
    value2_type: str
    value1_length: int
    value2_length: int
    row: int = location_field()
    index_column: str = location_field()
    index_value: str = location_field()
    position_run1: int = location_field()
    position_run2: int = location_field()
    column_index: int = location_field()
    column_name: Optional[str] = location_field()


@dataclass(frozen=True)
class HeaderFieldChangeDiscrepancy(DiscrepancyBase):
    """Represents a single header field change."""

    value1: Optional[str]
    value2: Optional[str]
    column_index: int = location_field()
    row: int = location_field()  # Headers are always row 1


@dataclass(frozen=True)
class ColumnReorderDiscrepancy(DiscrepancyBase):
    """Represents a single column being reordered."""

    column_name: str = location_field()
    position_run1: int = location_field()
    position_run2: int = location_field()


@dataclass(frozen=True)
class RowReorderDiscrepancy(DiscrepancyBase):
    """Represents a single row being reordered."""

    index_column: str = location_field()
    index_value: str = location_field()
    position_run1: int = location_field()
    position_run2: int = location_field()


# Union type for all specific discrepancy classes
Discrepancy = Union[
    MissingFileDiscrepancy,
    MissingDirectoryDiscrepancy,
    RowCountDifferenceDiscrepancy,
    ColumnCountDifferenceDiscrepancy,
    FileReadErrorDiscrepancy,
    EmptyFileDiscrepancy,
    NoIndexColumnDiscrepancy,
    DuplicateColumnNamesDiscrepancy,
    MissingRowDiscrepancy,
    ExtraRowDiscrepancy,
    # New flat discrepancy types (preferred)
    FieldChangeDiscrepancy,
    HeaderFieldChangeDiscrepancy,
    ColumnReorderDiscrepancy,
    RowReorderDiscrepancy,
    # Legacy aggregated types (kept for backward compatibility, but flattened in practice)
    HeaderDifferenceDiscrepancy,
    RowDifferenceDiscrepancy,
    ColumnOrderDifferenceDiscrepancy,
    RowOrderDifferenceDiscrepancy,
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
