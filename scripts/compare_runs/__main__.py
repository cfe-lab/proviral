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
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Any, TypeAlias

from .comparison_report import ComparisonReport
from .errors import (
    FileReadError,
    NoIndexColumn,
    MultipleIndexColumns,
)
from .discrepancy import (
    Severity,
    Confidence,
    _trim_value_for_display,
    Discrepancy as ImportedDiscrepancy,
    # Import all specific discrepancy classes
    DuplicateColumnNames,
    ColumnCountDifference,
    RowCountDifference,
    MissingRow,
    ExtraRow,
    MissingFile,
    MissingDirectory,
    # New flat discrepancy types
    FieldChange,
    HeaderFieldChange,
    HeaderDifference,
    ColumnReorder,
    RowReorder,
)
from .csv_utils import (
    get_file_metadata,
    validate_column_names,
    create_column_mapping,
    analyze_column_differences,
    compare_column_orders,
)

# from .index_discovery import (
#     discover_index_column,  # Replaced with regex-based approach
# )
from .row_comparison import (
    analyze_header_differences,
    analyze_row_differences,
)
from .file_operations import (
    find_versions,
    get_proviral_csv_files,
    read_csv_file,
    _build_index_row_mapping,
)


# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


Discrepancy: TypeAlias = ImportedDiscrepancy  # type: ignore[assignment]


def compare_csv_contents(
    report: ComparisonReport,
    file1: Path,
    file2: Path,
    version: str,
    filename: str,
    index_pattern: str,
) -> int:
    """
    Compare the contents of two CSV files and return number of discrepancies and errors.

    Args:
        file1: Path to first CSV file
        file2: Path to second CSV file
        version: Version being compared
        filename: Name of the file being compared

    Returns:
        int: Number of discrepancies and errors found during comparison.
    """

    content1 = read_csv_file(file1)
    content2 = read_csv_file(file2)
    initial = len(report.results) + len(report.errors)

    # Check if either file is empty or failed to read
    if not content1 and not content2:
        return 0  # Both empty, no discrepancies

    if not content1:
        report.add_error(
            FileReadError(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description="First file is empty or failed to read",
                file=filename,
                version=version,
                run="run1",
                file_path=str(file1),
            )
        )
        return 1

    if not content2:
        report.add_error(
            FileReadError(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description="Second file is empty or failed to read",
                file=filename,
                version=version,
                run="run2",
                file_path=str(file2),
            )
        )
        return 1

    # Get headers if available
    headers1 = content1[0] if content1 else []
    headers2 = content2[0] if content2 else []

    # Validate column names for duplicates
    dup_check1 = validate_column_names(content1)
    if dup_check1:
        report.add_discrepancy(
            DuplicateColumnNames(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"Duplicate column name '{dup_check1['duplicate_header']}' at positions {dup_check1['positions']} in run1",
                file=filename,
                version=version,
                run="run1",
                all_headers=dup_check1["all_headers"],
                duplicate_header=dup_check1["duplicate_header"],
                duplicate_column=dup_check1["duplicate_header"],
                positions=dup_check1["positions"],
            ),
        )

    dup_check2 = validate_column_names(content2)
    if dup_check2:
        report.add_discrepancy(
            DuplicateColumnNames(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"Duplicate column name '{dup_check2['duplicate_header']}' at positions {dup_check2['positions']} in run2",
                file=filename,
                version=version,
                run="run2",
                all_headers=dup_check2["all_headers"],
                duplicate_header=dup_check2["duplicate_header"],
                duplicate_column=dup_check2["duplicate_header"],
                positions=dup_check2["positions"],
            ),
        )

    # Check for column order differences (only if no duplicate names and same columns)
    order_diff = compare_column_orders(headers1, headers2)
    if order_diff:
        # Create individual discrepancies for each reordered column
        for reordered_col in order_diff["reordered_columns"]:
            report.add_discrepancy(
                ColumnReorder(
                    severity=Severity.MEDIUM,
                    confidence=Confidence.HIGH,
                    description=f"Column '{reordered_col['column']}' moved from position {reordered_col['position_run1']} to {reordered_col['position_run2']}",
                    file=filename,
                    version=version,
                    run="both",
                    column_name=reordered_col["column"],
                    position_run1=reordered_col["position_run1"],
                    position_run2=reordered_col["position_run2"],
                ),
            )

    # Create column mappings for name-based access (if no duplicates)
    column_map1 = create_column_mapping(headers1) if not dup_check1 else {}
    column_map2 = create_column_mapping(headers2) if not dup_check2 else {}

    # Resolve index column using regex pattern matching
    index_result = resolve_index_column_by_regex(
        content1, content2, filename, index_pattern
    )

    # Log index column resolution results and handle different outcomes
    if index_result["status"] == "success":
        index_column_name = index_result["column_name"]
        logger.debug(
            f"Found index column for {filename}: {index_column_name} (matched pattern: {index_pattern})"
        )
    elif index_result["status"] == "multiple_matches":
        # Report error for multiple matching columns
        report.add_error(
            MultipleIndexColumns(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"Multiple columns match pattern '{index_pattern}' for file {filename}: {', '.join(index_result['matching_columns'])}. Please use a more specific pattern.",
                file=filename,
                version=version,
                run="both",
                pattern=index_pattern,
                matching_columns=index_result["matching_columns"],
                matching_pairs=index_result["matching_pairs"],
            )
        )
        return len(report.results) + len(report.errors) - initial
    elif index_result["status"] == "error":
        # Report error for regex or other issues
        report.add_error(
            FileReadError(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"Error resolving index column: {index_result['error_message']}",
                file=filename,
                version=version,
                run="both",
                file_path=f"Pattern: {index_pattern}",
            )
        )
        return len(report.results) + len(report.errors) - initial
    else:  # status == "no_match"
        index_column_name = None
        logger.debug(
            f"No index column found for {filename} with pattern: {index_pattern}"
        )

    # Check if files have same number of rows
    if len(content1) != len(content2):
        report.add_discrepancy(
            RowCountDifference(
                severity=Severity.HIGH,
                confidence=Confidence.HIGH,
                description=f"Row count differs: {len(content1)} vs {len(content2)} (difference: {abs(len(content1) - len(content2))} rows)",
                file=filename,
                version=version,
                run="both",
                run1_rows=len(content1),
                run2_rows=len(content2),
                row_difference=abs(len(content1) - len(content2)),
            ),
        )

    # Always check headers first, regardless of comparison method
    if headers1 != headers2:
        # Identify specific header differences
        changed_headers = analyze_header_differences(headers1, headers2)

        # Create individual discrepancies for each header field change
        for removed_header in changed_headers["removed"]:
            report.add_discrepancy(
                HeaderFieldChange(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"Header '{removed_header['value']}' removed at position {removed_header['index']}",
                    file=filename,
                    version=version,
                    run="run1",  # Present in run1, missing in run2
                    column_index=removed_header["index"],
                    value1=removed_header["value"],
                    value2=None,
                    row=1,  # Headers are always row 1
                ),
            )

        for added_header in changed_headers["added"]:
            report.add_discrepancy(
                HeaderFieldChange(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"Header '{added_header['value']}' added at position {added_header['index']}",
                    file=filename,
                    version=version,
                    run="run2",  # Missing in run1, present in run2
                    column_index=added_header["index"],
                    value1=None,
                    value2=added_header["value"],
                    row=1,  # Headers are always row 1
                ),
            )

        for modified_header in changed_headers["modified"]:
            report.add_discrepancy(
                HeaderFieldChange(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"Header changed at position {modified_header['index']}: '{modified_header['old_header']}' → '{modified_header['new_header']}'",
                    file=filename,
                    version=version,
                    run="both",
                    value1=modified_header["old_header"],
                    value2=modified_header["new_header"],
                    column_index=modified_header["index"],
                    row=1,  # Headers are always row 1
                ),
            )

        # Check for column count differences in header
        if len(headers1) != len(headers2):
            missing_cols, extra_cols = analyze_column_differences(headers1, headers2)
            report.add_discrepancy(
                HeaderDifference(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"Header column count differs: {len(headers1)} vs {len(headers2)} ({len(missing_cols)} missing, {len(extra_cols)} extra in run2)",
                    file=filename,
                    version=version,
                    run="both",
                    columns_missing_in_run2=missing_cols,
                    columns_extra_in_run2=extra_cols,
                    run1_header_count=len(headers1),
                    run2_header_count=len(headers2),
                ),
            )

    # Use index-column-based comparison if available, otherwise report error
    if index_column_name is not None and not dup_check1 and not dup_check2:
        # Use index-column-based row comparison (skip header row since already compared)
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
        for discrepancy in index_based_discrepancies:
            report.add_discrepancy(discrepancy)
    else:
        # Report error when no index column is available - cannot perform row comparison
        logger.warning(
            f"Cannot perform row comparison for {filename}: no suitable index column available"
        )

        # Determine the specific reason for unavailability
        reason = "unknown"
        if dup_check1 or dup_check2:
            reason = "duplicate_column_names"
        elif not index_column_name:
            reason = "no_regex_match"
        else:
            reason = "no_index_available"

        report.add_error(
            NoIndexColumn(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"Row comparison skipped: No suitable index column available ({reason})",
                file=filename,
                version=version,
                run="both",
                reason=reason,
                index_analysis={
                    "pattern": index_pattern,
                    "column_found": index_column_name,
                },
                has_duplicate_columns_run1=bool(dup_check1),
                has_duplicate_columns_run2=bool(dup_check2),
            )
        )

    return len(report.results) + len(report.errors) - initial


def _determine_field_change_severity(change_type: str) -> Severity:
    """Determine severity of individual field change based on change type."""
    if change_type == "outcome":
        return Severity.HIGH
    elif change_type == "numeric":
        return Severity.MEDIUM
    elif change_type == "numeric_equivalent":
        return Severity.LOW
    else:
        return Severity.MEDIUM  # Default for text changes


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

    # Check for numeric equivalent changes (lower severity than actual numeric changes)
    numeric_equivalent_changes = sum(
        1
        for change in column_differences.get("field_changes", [])
        if change["change_type"] == "numeric_equivalent"
    )
    if numeric_equivalent_changes > 0:
        return Severity.LOW

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
    discrepancies: List[Discrepancy] = []

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
                column_differences = analyze_row_differences(
                    row1, row2, headers1, headers2, column_map1, column_map2
                )

                # Create individual discrepancies for each field change
                for field_change in column_differences["field_changes"]:
                    field_severity = _determine_field_change_severity(
                        field_change["change_type"]
                    )
                    field_confidence = _determine_row_difference_confidence(row1, row2)

                    col_name = (
                        field_change["column_name"]
                        or f"column_{field_change['column_index']}"
                    )

                    discrepancies.append(
                        FieldChange(
                            severity=field_severity,
                            confidence=field_confidence,
                            description=f"Field '{col_name}' differs in row with {index_column_name}='{_trim_value_for_display(str(index_value))}': '{_trim_value_for_display(field_change['value1'])}' → '{_trim_value_for_display(field_change['value2'])}'",
                            file=filename,
                            version=version,
                            run="both",
                            row=pos1 + 1,  # 1-based row number
                            index_column=index_column_name,
                            index_value=str(index_value),
                            position_run1=pos1,
                            position_run2=pos2,
                            column_index=field_change["column_index"],
                            column_name=field_change["column_name"],
                            change_type=field_change["change_type"],
                            value1=field_change["value1"],
                            value2=field_change["value2"],
                            value1_type=field_change["value1_type"],
                            value2_type=field_change["value2_type"],
                            value1_length=field_change["value1_length"],
                            value2_length=field_change["value2_length"],
                        )
                    )

                # Check for column count differences in data rows
                if len(row1) != len(row2):
                    # For data rows, missing_cols and extra_cols will be []
                    missing_cols, extra_cols = analyze_column_differences(row1, row2)
                    discrepancies.append(
                        ColumnCountDifference(
                            severity=Severity.HIGH,
                            confidence=Confidence.HIGH,
                            description=f"Row with {index_column_name}='{_trim_value_for_display(str(index_value))}' column count differs: {len(row1)} vs {len(row2)} ({len(missing_cols)} missing, {len(extra_cols)} extra in run2)",
                            file=filename,
                            version=version,
                            run="both",
                            run1_columns=len(row1),
                            run2_columns=len(row2),
                            columns_missing_in_run2=missing_cols,
                            columns_extra_in_run2=extra_cols,
                            row=pos1 + 1,  # 1-based row number
                            column_difference=abs(len(row1) - len(row2)),
                        )
                    )

        elif row1_info and not row2_info:
            # Row exists only in first file
            row1 = row1_info["row"]
            pos1 = row1_info["position"]

            discrepancies.append(
                MissingRow(
                    severity=Severity.HIGH,
                    confidence=Confidence.HIGH,
                    description=f"Row with {index_column_name}='{_trim_value_for_display(str(index_value))}' missing in run2",
                    file=filename,
                    version=version,
                    run="run2",
                    index_column=index_column_name,
                    index_value=str(index_value),
                    position_run1=pos1,
                    missing_from="run2",
                    present_in="run1",
                    missing_row_data=row1,
                )
            )

        elif row2_info and not row1_info:
            # Row exists only in second file
            row2 = row2_info["row"]
            pos2 = row2_info["position"]

            discrepancies.append(
                ExtraRow(
                    severity=Severity.HIGH,
                    confidence=Confidence.HIGH,
                    description=f"Row with {index_column_name}='{_trim_value_for_display(str(index_value))}' extra in run2",
                    file=filename,
                    version=version,
                    run="run2",
                    index_column=index_column_name,
                    index_value=str(index_value),
                    position_run2=pos2,
                    missing_from="run1",
                    present_in="run2",
                    extra_row_data=row2,
                )
            )

    # Report row order differences if rows exist in both files but in different positions
    if position_differences:
        # Create individual discrepancies for each reordered row
        for position_diff in position_differences:
            discrepancies.append(
                RowReorder(
                    severity=Severity.LOW,
                    confidence=Confidence.HIGH,
                    description=f"Row with {index_column_name}='{_trim_value_for_display(position_diff['index_value'])}' moved from position {position_diff['position_run1']} to {position_diff['position_run2']}",
                    file=filename,
                    version=version,
                    run="both",
                    index_column=index_column_name,
                    index_value=position_diff["index_value"],
                    position_run1=position_diff["position_run1"],
                    position_run2=position_diff["position_run2"],
                )
            )

    return discrepancies


def compare_proviral_files(
    run1_dir: Path,
    run2_dir: Path,
    version: str,
    report: ComparisonReport,
    index_pattern: str,
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
            MissingDirectory(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description=f"No proviral CSV files found for {version}",
                file="proviral_directory",
                version=version,
                run="both",
                missing_from="both",
                present_in="",
            ),
        )
        return

    for csv_name in sorted(all_csv_names):
        if csv_name not in csv_files1:
            # Get file info from run2
            file2_info = get_file_metadata(csv_files2[csv_name])
            report.add_discrepancy(
                MissingFile(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"File missing in run1: {csv_name} (run2 size: {file2_info['size']} bytes)",
                    file=csv_name,
                    version=version,
                    run="run1",
                    missing_from="run1",
                    present_in="run2",
                ),
            )
            continue

        if csv_name not in csv_files2:
            # Get file info from run1
            file1_info = get_file_metadata(csv_files1[csv_name])
            report.add_discrepancy(
                MissingFile(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.HIGH,
                    description=f"File missing in run2: {csv_name} (run1 size: {file1_info['size']} bytes)",
                    file=csv_name,
                    version=version,
                    run="run2",
                    missing_from="run2",
                    present_in="run1",
                ),
            )
            continue

        # Compare the files
        count = compare_csv_contents(
            report,
            csv_files1[csv_name],
            csv_files2[csv_name],
            version,
            csv_name,
            index_pattern,
        )

        if not count:
            report.mark_file_identical()


def compare_runs(
    run1_dir: Path, run2_dir: Path, index_pattern: str
) -> ComparisonReport:
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
            MissingDirectory(
                severity=Severity.CRITICAL,
                confidence=Confidence.HIGH,
                description="No common versions found between the two runs",
                file="versions",
                version="all",
                run="both",
                missing_from="both",
                present_in="",
            ),
        )
        return report

    # Report version mismatches
    only_in_run1 = set(versions1) - set(versions2)
    only_in_run2 = set(versions2) - set(versions1)

    if only_in_run1:
        report.add_discrepancy(
            MissingDirectory(
                severity=Severity.MEDIUM,
                confidence=Confidence.HIGH,
                description=f"Versions only in run1: {sorted(only_in_run1)}",
                file="versions",
                version="multiple",
                run="run1",
                missing_from="run2",
                present_in="run1",
            ),
        )

    if only_in_run2:
        report.add_discrepancy(
            MissingDirectory(
                severity=Severity.MEDIUM,
                confidence=Confidence.HIGH,
                description=f"Versions only in run2: {sorted(only_in_run2)}",
                file="versions",
                version="multiple",
                run="run2",
                missing_from="run1",
                present_in="run2",
            ),
        )

    # Compare each common version
    for version in sorted(common_versions):
        try:
            compare_proviral_files(run1_dir, run2_dir, version, report, index_pattern)
        except Exception as e:
            report.add_error(
                FileReadError(
                    severity=Severity.CRITICAL,
                    confidence=Confidence.MEDIUM,
                    description=f"Error comparing {version}: {str(e)}",
                    file="comparison_error",
                    version=version,
                    run="both",
                    file_path=str(e),
                ),
            )

    return report


def resolve_index_column_by_regex(
    content1: List[List[str]],
    content2: List[List[str]],
    filename: str,
    index_pattern: str,
) -> Dict[str, Any]:
    """
    Resolve index column using regex pattern matching against filename/columnname.

    Args:
        content1: First CSV data
        content2: Second CSV data
        filename: Name of the CSV file being compared
        index_pattern: Regex pattern to match against filename/columnname

    Returns:
        Dictionary with resolution results:
        - status: "success", "no_match", "multiple_matches", "error"
        - column_name: str (if status == "success")
        - matching_columns: List[str] (if status == "multiple_matches")
        - matching_pairs: List[str] (if status == "multiple_matches")
        - error_message: str (if status == "error")
    """
    if not content1 or not content2:
        return {"status": "error", "error_message": "Empty or missing CSV data"}

    headers1 = content1[0] if content1 else []
    headers2 = content2[0] if content2 else []

    # Compile the regex pattern
    try:
        pattern = re.compile(index_pattern)
    except re.error as e:
        logger.error(f"Invalid regex pattern: {index_pattern}")
        return {"status": "error", "error_message": f"Invalid regex pattern: {e}"}

    # Check each column name against the pattern and collect all matches
    # We need to find columns that exist in both files and match the pattern
    matching_columns = []
    matching_pairs = []

    # Find columns that exist in both files
    common_columns = set(headers1) & set(headers2)

    for col_name in common_columns:
        match_string = f"{filename}/{col_name}"
        if pattern.search(match_string):
            matching_columns.append(col_name)
            matching_pairs.append(match_string)
            logger.debug(
                f"Index column match found: {match_string} matches pattern {index_pattern}"
            )

    if len(matching_columns) == 0:
        logger.debug(
            f"No index column match found for file {filename} with pattern {index_pattern}"
        )
        return {"status": "no_match"}
    elif len(matching_columns) == 1:
        return {
            "status": "success",
            "column_name": matching_columns[0],
            "match_string": matching_pairs[0],
        }
    else:
        logger.warning(
            f"Multiple index columns match pattern {index_pattern} for file {filename}: {matching_columns}"
        )
        return {
            "status": "multiple_matches",
            "matching_columns": matching_columns,
            "matching_pairs": matching_pairs,
        }


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Compare proviral analysis results between different pipeline runs",
        epilog="Example: ws-compare-runs --indexes '.*/sample' /path/to/run1 /path/to/run2",
    )

    parser.add_argument(
        "--indexes",
        default=".*/sample|.*/samp_name",
        help="Regex pattern to match index columns in format 'filename/columnname'. Examples: '.*/sample' or 'file1.csv/id|file2.csv/sample_id'",
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
        report = compare_runs(args.run1_dir, args.run2_dir, args.indexes)

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
