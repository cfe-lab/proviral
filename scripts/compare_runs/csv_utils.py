"""
CSV utility functions for the compare_runs script.

This module provides utility functions for working with CSV files,
including validation, metadata extraction, and data type detection.
"""

from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any


def get_file_metadata(file_path: Path) -> Dict[str, Any]:
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


def validate_column_names(csv_data: List[List[str]]) -> Optional[Dict[str, Any]]:
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


def create_column_mapping(headers: List[str]) -> Dict[str, int]:
    """
    Create a mapping from column names to their indices.

    Args:
        headers: List of column names

    Returns:
        Dictionary mapping column names to indices
    """
    return {header: i for i, header in enumerate(headers)}


def get_column_value_by_name(
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


def is_numeric(value: str) -> bool:
    """Check if a string represents a numeric value."""
    try:
        float(value)
        return True
    except ValueError:
        return False


def get_value_type(value: str) -> str:
    """Get the type classification of a value."""
    if not value:
        return "empty"
    elif is_numeric(value):
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


def classify_value_change(val1: str, val2: str) -> str:
    """Classify the type of change between two values."""
    if not val1 and val2:
        return "added"
    elif val1 and not val2:
        return "removed"
    elif is_numeric(val1) and is_numeric(val2):
        # Check if numeric values are equal despite string differences
        try:
            if float(val1) == float(val2):
                return "numeric_equivalent"  # Same numeric value, different string representation
            else:
                return "numeric"  # Different numeric values
        except ValueError:
            # Fallback if conversion fails (shouldn't happen since is_numeric passed)
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


def analyze_column_differences(row1: List[str], row2: List[str]) -> tuple[int, int]:
    """Analyze column count differences between two rows."""
    len1, len2 = len(row1), len(row2)
    if len1 > len2:
        return len1 - len2, 0  # missing in run2, extra in run1
    elif len2 > len1:
        return 0, len2 - len1  # missing in run1, extra in run2
    else:
        return 0, 0


def compare_column_orders(
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
