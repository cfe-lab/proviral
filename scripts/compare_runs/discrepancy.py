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

# Import core functionality
from dataclasses import dataclass, fields
from typing import Any, Dict, List, Optional, Type, Union, TypeAlias
from .severity import Severity
from .confidence import Confidence


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


# Helper function to mark trimmable fields on dataclasses
def _mark_trimmable_fields(cls: Type, trimmable: List[str]) -> Type:
    """
    Mark trimmable fields in the class metadata.

    Args:
        cls: The dataclass to mark
        trimmable: List of field names that should be trimmed for display
    """

    # Mark trimmable fields in the class metadata
    if not hasattr(cls, "_trimmable_fields"):
        cls._trimmable_fields = set(trimmable)
    else:
        cls._trimmable_fields.update(trimmable)

    return cls


def location_fields(locations: List[str]):
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


def trimmable_fields(trimmable: List[str]):
    """
    Decorator to mark trimmable fields on a dataclass.

    This decorator should be applied before @dataclass to mark which fields
    should be trimmed for display purposes when serializing to dict.

    Args:
        trimmable: List of field names that should be trimmed for display

    Example:
        @dataclass(frozen=True)
        class MyDiscrepancy(DiscrepancyBase):
            value1: str
            value2: str
            other_field: str
    """

    def decorator(cls: Type) -> Type:
        return _mark_trimmable_fields(cls, trimmable)

    return decorator


def _trim_value_for_display(value: str, max_length: int = 50) -> str:
    """
    Trim a value for display purposes if it exceeds the maximum length.

    For values longer than max_length, truncates with "..." in the middle
    to preserve both the beginning and end of the value.

    Args:
        value: The string value to potentially trim
        max_length: Maximum allowed length (default: 50)

    Returns:
        The original value if <= max_length, otherwise a trimmed version

    Examples:
        >>> _trim_value_for_display("short")
        'short'
        >>> _trim_value_for_display("this_is_a_very_long_value_that_exceeds_fifty_characters")
        'this_is_a_...characters'
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


Serialized: TypeAlias = Union[int, float, bool, str, List["Serialized"]]


def serialize_value(value: object) -> Union[Serialized]:
    """
    Serialize a value for JSON output.
    Converts basic types to their JSON-compatible representations.
    """

    if isinstance(value, (int, float, bool)):
        return value
    elif isinstance(value, list):
        return [serialize_value(item) for item in value]
    else:
        return str(value)


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
        top_level_dict: Dict[str, Any] = {
            "type": self.__class__.__name__,
        }

        # Get location fields and trimmable fields from class metadata
        location_fields: set = getattr(self.__class__, "_location_fields", set())
        trimmable_fields: set = getattr(self.__class__, "_trimmable_fields", set())

        # Traverse all dataclass fields
        for field_info in fields(self):
            value = getattr(self, field_info.name)

            # Apply trimming if this field is marked as trimmable
            if field_info.name in trimmable_fields and isinstance(value, str):
                value = _trim_value_for_display(value)

            # Check if this is a location field
            if field_info.name in location_fields:
                # Add to location dict
                location_dict[field_info.name] = serialize_value(value)
            else:
                # Add to top level
                top_level_dict[field_info.name] = serialize_value(value)

        top_level_dict["location"] = location_dict

        return top_level_dict


@dataclass(frozen=True)
class MissingFile(DiscrepancyBase):
    """Represents a missing file discrepancy."""

    missing_from: str
    present_in: str


@dataclass(frozen=True)
class MissingDirectory(DiscrepancyBase):
    """Represents a missing directory discrepancy."""

    missing_from: str
    present_in: str


@location_fields(["row", "changed_headers", "total_header_changes"])
@dataclass(frozen=True)
class HeaderDifference(DiscrepancyBase):
    """Represents a header difference discrepancy."""

    header_changes: Dict[str, Any]
    run1_header_count: int
    run2_header_count: int
    row: int
    changed_headers: List[str]
    total_header_changes: int


@trimmable_fields(["index_value"])
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
class RowDifference(DiscrepancyBase):
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
class RowCountDifference(DiscrepancyBase):
    """Represents a row count difference discrepancy."""

    run1_rows: int
    run2_rows: int
    row_difference: int


@location_fields(
    [
        "row",
    ]
)
@dataclass(frozen=True)
class ColumnCountDifference(DiscrepancyBase):
    """Represents a column count difference discrepancy."""

    run1_columns: int
    run2_columns: int
    columns_missing_in_run2: List[str]
    columns_extra_in_run2: List[str]
    row: int
    column_difference: int


@location_fields(["positions"])
@dataclass(frozen=True)
class DuplicateColumnNames(DiscrepancyBase):
    """Represents a duplicate column names discrepancy."""

    all_headers: List[str]
    duplicate_header: str  # alias for duplicate_column
    duplicate_column: str  # Original field name
    positions: List[int]


@location_fields(["reordered_columns"])
@dataclass(frozen=True)
class ColumnOrderDifference(DiscrepancyBase):
    """Represents a column order difference discrepancy."""

    order_differences: Dict[str, Any]
    reordered_count: int
    reordered_columns: List[str]


@location_fields(["index_column", "reordered_rows"])
@dataclass(frozen=True)
class RowOrderDifference(DiscrepancyBase):
    """Represents a row order difference discrepancy."""

    position_differences: List[Dict[str, Any]]
    reordered_count: int
    index_column: str
    reordered_rows: List[str]


@trimmable_fields(["index_value"])
@location_fields(
    ["index_column", "index_value", "position_run1", "missing_from", "present_in"]
)
@dataclass(frozen=True)
class MissingRow(DiscrepancyBase):
    """Represents a missing row discrepancy."""

    missing_row_data: List[str]
    index_column: str
    index_value: str
    position_run1: int
    missing_from: str
    present_in: str


@trimmable_fields(["index_value"])
@location_fields(
    ["index_column", "index_value", "position_run1", "missing_from", "present_in"]
)
@dataclass(frozen=True)
class ExtraRow(DiscrepancyBase):
    """Represents an extra row discrepancy."""

    extra_row_data: List[str]
    index_column: str
    index_value: str
    position_run2: int
    missing_from: str
    present_in: str


@trimmable_fields(["value1", "value2", "index_value"])
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
class FieldChange(DiscrepancyBase):
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
class HeaderFieldChange(DiscrepancyBase):
    """Represents a single header field change."""

    value1: Optional[str]
    value2: Optional[str]
    column_index: int
    row: int  # Headers are always row 1


@location_fields(["column_name", "position_run1", "position_run2"])
@dataclass(frozen=True)
class ColumnReorder(DiscrepancyBase):
    """Represents a single column being reordered."""

    column_name: str
    position_run1: int
    position_run2: int


@location_fields(["index_column", "index_value", "position_run1", "position_run2"])
@dataclass(frozen=True)
class RowReorder(DiscrepancyBase):
    """Represents a single row being reordered."""

    index_column: str
    index_value: str
    position_run1: int
    position_run2: int


# Union type for all specific discrepancy classes
Discrepancy = Union[
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
]
