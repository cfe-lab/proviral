"""
Discrepancy classes and enums for the compare_runs script.

This module defines the core classes and enums used to represent and manage
discrepancies found during proviral analysis comparison.
"""

import json
from collections import defaultdict
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional


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
            "description": self.description,  # Descriptions are usually crafted and might not need aggressive trimming here
            "location": _trim_data_recursively(
                self.location
            ),  # Trim location data as well
            "values": _trim_data_recursively(self.values),  # Recursively trim values
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
