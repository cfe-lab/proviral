"""
Index column discovery functionality for CSV comparison.

This module provides utilities to automatically discover suitable index columns
for row-by-row comparison of CSV files, enabling more accurate comparison when
row order differs between files.
"""

from typing import List, Dict, Optional, Any


def find_unique_value_columns(csv_data: List[List[str]]) -> List[str]:
    """
    Find columns that contain only unique values (no duplicates).

    Args:
        csv_data: List of rows, where first row contains headers

    Returns:
        List of column names that have unique values
    """
    if not csv_data or len(csv_data) < 2:  # Need at least header + 1 data row
        return []

    headers = csv_data[0]
    unique_columns = []

    for col_idx, header in enumerate(headers):
        # Collect all non-empty values in this column (skip header row)
        values = []
        for row_idx in range(1, len(csv_data)):
            if col_idx < len(csv_data[row_idx]):
                value = csv_data[row_idx][col_idx].strip()
                if value:  # Only consider non-empty values
                    values.append(value)

        # Check if all values are unique
        if len(values) == len(set(values)):
            unique_columns.append(header)

    return unique_columns


def _count_shared_values(
    csv_data1: List[List[str]], csv_data2: List[List[str]], column_name: str
) -> int:
    """
    Count how many values are shared between the same column in two CSV files.

    Args:
        csv_data1: First CSV data
        csv_data2: Second CSV data
        column_name: Name of column to compare

    Returns:
        Number of shared non-empty values
    """
    if not csv_data1 or not csv_data2 or len(csv_data1) < 2 or len(csv_data2) < 2:
        return 0

    headers1 = csv_data1[0]
    headers2 = csv_data2[0]

    if column_name not in headers1 or column_name not in headers2:
        return 0

    col_idx1 = headers1.index(column_name)
    col_idx2 = headers2.index(column_name)

    # Collect non-empty values from both files
    values1 = set()
    for row_idx in range(1, len(csv_data1)):
        if col_idx1 < len(csv_data1[row_idx]):
            value = csv_data1[row_idx][col_idx1].strip()
            if value:
                values1.add(value)

    values2 = set()
    for row_idx in range(1, len(csv_data2)):
        if col_idx2 < len(csv_data2[row_idx]):
            value = csv_data2[row_idx][col_idx2].strip()
            if value:
                values2.add(value)

    # Count shared values
    return len(values1 & values2)


def count_shared_values(
    csv_data1: List[List[str]], csv_data2: List[List[str]], column_name: str
) -> int:
    """
    Count how many values are shared between the same column in two CSV files.

    This is the public interface to _count_shared_values.

    Args:
        csv_data1: First CSV data
        csv_data2: Second CSV data
        column_name: Name of column to compare

    Returns:
        Number of shared non-empty values
    """
    return _count_shared_values(csv_data1, csv_data2, column_name)


def discover_index_column(
    csv_data1: List[List[str]], csv_data2: List[List[str]]
) -> Optional[Dict[str, Any]]:
    """
    Discover the best column to use as an index for row comparison.

    Args:
        csv_data1: First CSV data
        csv_data2: Second CSV data

    Returns:
        Dictionary with index column information or None if no suitable column found.

        Success case returns:
        {
            "index_column": int,        # Column index (0-based)
            "column_name": str,         # Column name
            "shared_values": int,       # Number of shared values
            "reason": "success",
            "column_analysis": {...}    # Detailed analysis
        }

        Failure cases return:
        - None: No unique columns in both files
        - Dict with index_column=None: Unique columns exist but no shared values
    """
    if not csv_data1 or not csv_data2 or len(csv_data1) < 2 or len(csv_data2) < 2:
        return None

    headers1 = csv_data1[0]
    headers2 = csv_data2[0]

    # Headers must match for safe column comparison
    if headers1 != headers2:
        return None

    # Find columns that are unique in both files
    unique_cols1 = set(find_unique_value_columns(csv_data1))
    unique_cols2 = set(find_unique_value_columns(csv_data2))
    common_unique_cols = unique_cols1 & unique_cols2

    if not common_unique_cols:
        # No columns are unique in both files
        return None

    # Evaluate each candidate column by counting shared values
    best_column = None
    best_shared_count = 0
    best_column_name = None
    column_analysis = {}

    # Sort column names for consistent tie-breaking
    for col_name in sorted(common_unique_cols):
        shared_count = _count_shared_values(csv_data1, csv_data2, col_name)
        column_analysis[col_name] = {
            "shared_values": shared_count,
            "column_index": headers1.index(col_name),
        }

        if shared_count > best_shared_count:
            best_shared_count = shared_count
            best_column = headers1.index(col_name)
            best_column_name = col_name

    if best_shared_count == 0:
        # Unique columns exist but no shared values
        return {
            "index_column": None,
            "shared_values": 0,
            "reason": "no_shared_values",
            "candidate_columns": len(common_unique_cols),
            "column_analysis": column_analysis,
        }

    # Success case
    return {
        "index_column": best_column,
        "column_name": best_column_name,
        "shared_values": best_shared_count,
        "reason": "success",
        "column_analysis": column_analysis,
    }
