"""
File operations module for compare_runs script.

This module contains functions for file system operations including:
- Finding version directories
- Getting CSV files from proviral subfolders
- Reading CSV file contents
- Building index row mappings
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Any

logger = logging.getLogger(__name__)


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


def _build_index_row_mapping(
    csv_data: List[List[str]], index_column_name: str
) -> Dict[str, Dict[str, Any]]:
    """
    Build a mapping from index column values to row data.

    Args:
        csv_data: CSV data as list of lists
        index_column_name: Name of column to use as index

    Returns:
        Dictionary mapping index values to row data and position
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
