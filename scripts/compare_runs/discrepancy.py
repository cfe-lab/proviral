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

from dataclasses import dataclass, fields
from typing import Any, Dict, List, Optional, Union, Type
from enum import Enum


# Helper function to mark location fields on dataclasses
def _mark_location_fields(cls: Type, locations: Optional[List[str]] = None) -> Type:
    """
    Mark location fields in the class metadata.

    Args:
        cls: The dataclass to mark
        locations: List of field names that should be treated as location fields
    """
    if locations is None:
        locations = []

    # Mark location fields in the class metadata
    if not hasattr(cls, "_location_fields"):
        cls._location_fields = set(locations)
    else:
        cls._location_fields.update(locations)

    return cls


def location_fields(locations: Optional[List[str]] = None):
    """
    Decorator to mark location fields on a dataclass.

    This decorator should be applied before @dataclass to mark which fields
    should be treated as location fields in the serialized output.

    Args:
        locations: List of field names that should be treated as location fields

    Example:
        @location_fields(["row", "column_index"])
        @dataclass(frozen=True)
        class MyDiscrepancy(DiscrepancyBase):
            row: int
            column_index: int
            other_field: str
    """

    def decorator(cls: Type) -> Type:
        return _mark_location_fields(cls, locations)

    return decorator


# Severity levels for discrepancies
class Severity(Enum):
    CRITICAL = "CRITICAL"
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"

    def __str__(self) -> str:
        return self.value


# Confidence levels for discrepancies
class Confidence(Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"

    def __str__(self) -> str:
        return self.value


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
    elif isinstance(data, Enum):
        return data.value  # Use the enum value directly
    elif isinstance(data, dict):
        return {
            key: _trim_data_recursively(value, max_length)
            for key, value in data.items()
        }
    return data


# Base class for all discrepancy types
@location_fields(["file", "version", "run"])
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

        # Get location fields from class metadata
        location_fields: set = getattr(self.__class__, "_location_fields", set())

        # Traverse all dataclass fields
        for field_info in fields(self):
            value = getattr(self, field_info.name)

            # Skip None values for optional fields
            if value is None:
                continue

            # Check if this is a location field
            if field_info.name in location_fields:
                # Add to location dict
                location_dict[field_info.name] = value
            else:
                # Add to top level
                top_level_dict[field_info.name] = value

        top_level_dict["location"] = location_dict
        trimmed = _trim_data_recursively(top_level_dict)
        return {
            **trimmed,
            "type": self.__class__.__name__,
            "description": self.description,
        }


@location_fields()
@dataclass(frozen=True)
class MissingFileDiscrepancy(DiscrepancyBase):
    """Represents a missing file discrepancy."""

    missing_from: str
    present_in: str


@location_fields()
@dataclass(frozen=True)
class MissingDirectoryDiscrepancy(DiscrepancyBase):
    """Represents a missing directory discrepancy."""

    missing_from: str
    present_in: str


@location_fields(["row", "changed_headers", "total_header_changes"])
@dataclass(frozen=True)
class HeaderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a header difference discrepancy."""

    header_changes: Dict[str, Any]
    run1_header_count: int
    run2_header_count: int
    row: int
    changed_headers: List[str]
    total_header_changes: int


@location_fields(
    [
        "row",
        "index_column",
        "index_value",
        "position_run1",
        "position_run2",
        "changed_columns",
        "total_field_changes",
    ]
)
@dataclass(frozen=True)
class RowDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row difference discrepancy."""

    field_changes: Dict[str, Any]
    change_types: List[str]
    row: int
    index_column: str
    index_value: str
    position_run1: int
    position_run2: int
    changed_columns: List[str]
    total_field_changes: int


@location_fields(["row_difference"])
@dataclass(frozen=True)
class RowCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row count difference discrepancy."""

    run1_rows: int
    run2_rows: int
    row_difference: int


@location_fields(
    [
        "row",
        "index_column",
        "index_value",
        "position_run1",
        "position_run2",
        "column_difference",
    ]
)
@dataclass(frozen=True)
class ColumnCountDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column count difference discrepancy."""

    run1_columns: int
    run2_columns: int
    columns_missing_in_run2: int
    columns_extra_in_run2: int
    row: int
    column_difference: int


@location_fields(["file_path"])
@dataclass(frozen=True)
class FileReadErrorDiscrepancy(DiscrepancyBase):
    """Represents a file read error discrepancy."""

    file_path: str


@location_fields(["file_path"])
@dataclass(frozen=True)
class EmptyFileDiscrepancy(DiscrepancyBase):
    """Represents an empty file discrepancy."""

    file_path: str
    file_size: int


@dataclass(frozen=True)
class NoIndexColumnDiscrepancy(DiscrepancyBase):
    """Represents a no index column discrepancy."""

    reason: str
    index_analysis: Optional[Dict[str, Any]]
    has_duplicate_columns_run1: bool
    has_duplicate_columns_run2: bool


@location_fields(["positions"])
@dataclass(frozen=True)
class DuplicateColumnNamesDiscrepancy(DiscrepancyBase):
    """Represents a duplicate column names discrepancy."""

    all_headers: List[str]
    duplicate_header: str  # alias for duplicate_column
    duplicate_column: str  # Original field name
    positions: List[int]


@location_fields(["reordered_columns"])
@dataclass(frozen=True)
class ColumnOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a column order difference discrepancy."""

    order_differences: Dict[str, Any]
    reordered_count: int
    reordered_columns: List[str]


@location_fields(["index_column", "reordered_rows"])
@dataclass(frozen=True)
class RowOrderDifferenceDiscrepancy(DiscrepancyBase):
    """Represents a row order difference discrepancy."""

    position_differences: List[Dict[str, Any]]
    reordered_count: int
    index_column: str
    reordered_rows: List[str]


@location_fields(
    ["index_column", "index_value", "position_run1", "missing_from", "present_in"]
)
@dataclass(frozen=True)
class MissingRowDiscrepancy(DiscrepancyBase):
    """Represents a missing row discrepancy."""

    missing_row_data: List[str]
    index_column: str
    index_value: str
    position_run1: int
    missing_from: str
    present_in: str


@location_fields(
    ["index_column", "index_value", "position_run1", "missing_from", "present_in"]
)
@dataclass(frozen=True)
class ExtraRowDiscrepancy(DiscrepancyBase):
    """Represents an extra row discrepancy."""

    extra_row_data: List[str]
    index_column: str
    index_value: str
    position_run2: int
    missing_from: str
    present_in: str


@location_fields(
    [
        "row",
        "index_column",
        "index_value",
        "position_run1",
        "position_run2",
        "column_index",
        "column_name",
    ],
)
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
    row: int
    index_column: str
    index_value: str
    position_run1: int
    position_run2: int
    column_index: int
    column_name: Optional[str]


@location_fields(["row", "column_index"])
@dataclass(frozen=True)
class HeaderFieldChangeDiscrepancy(DiscrepancyBase):
    """Represents a single header field change."""

    value1: Optional[str]
    value2: Optional[str]
    column_index: int
    row: int  # Headers are always row 1


@location_fields(["column_name", "position_run1", "position_run2"])
@dataclass(frozen=True)
class ColumnReorderDiscrepancy(DiscrepancyBase):
    """Represents a single column being reordered."""

    column_name: str
    position_run1: int
    position_run2: int


@location_fields(["index_column", "index_value", "position_run1", "position_run2"])
@dataclass(frozen=True)
class RowReorderDiscrepancy(DiscrepancyBase):
    """Represents a single row being reordered."""

    index_column: str
    index_value: str
    position_run1: int
    position_run2: int


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
