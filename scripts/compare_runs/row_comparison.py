"""
Row comparison functionality for CSV comparison.

This module provides utilities to analyze and summarize differences between
CSV rows and headers, supporting detailed comparison reporting.
"""

from typing import Dict, List, Optional, Any
from .discrepancy import _trim_value_for_display
from .csv_utils import classify_value_change, get_value_type


def get_header_change_summary(changes: Dict[str, Any]) -> str:
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


def analyze_row_differences(
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
                change_type = classify_value_change(val1, val2)
                changes["change_types"].append(change_type)

                change_info = {
                    "column_index": idx1,
                    "column_name": col_name,
                    "change_type": change_type,
                    "value1": _trim_value_for_display(val1),
                    "value2": _trim_value_for_display(val2),
                    "value1_type": get_value_type(val1),
                    "value2_type": get_value_type(val2),
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
                change_type = classify_value_change(val1, val2)
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
                    "value1": _trim_value_for_display(val1),
                    "value2": _trim_value_for_display(val2),
                    "value1_type": get_value_type(val1),
                    "value2_type": get_value_type(val2),
                    "value1_length": len(val1),
                    "value2_length": len(val2),
                }

                changes["field_changes"].append(change_info)

    return changes


def extract_column_names_from_changes(changes: Dict[str, Any]) -> List[str]:
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


def extract_header_names_from_changes(
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


def get_row_change_summary(changes: Dict[str, Any]) -> str:
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
