import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional


class ComparisonReport:
    """Collects and manages all discrepancies found during comparison."""

    def __init__(self, run1_dir: Path, run2_dir: Path):
        self.run1_dir = run1_dir
        self.run2_dir = run2_dir
        self.timestamp = datetime.now().isoformat()
        self.versions_run1: List[str] = []
        self.versions_run2: List[str] = []
        self.common_versions: List[str] = []
        self.results: List[Dict[str, Any]] = []
        # Processing errors during comparison (exceptions, etc.)
        self.errors: List[Dict[str, Any]] = []

    def add_discrepancy(self, version: str, filename: str, discrepancy):
        """Add a discrepancy to the report."""
        self.results.append(discrepancy.to_dict())

    def add_error(self, error):
        """Add a processing error to the report."""
        self.errors.append(error.to_dict())

    def mark_file_identical(
        self, version: str, filename: str, index_info: Optional[Dict[str, Any]] = None
    ):
        """Mark a file as identical between runs."""
        # In the flat structure, we don't track identical files separately
        # since they don't generate discrepancies
        pass

    def set_file_index_info(
        self, version: str, filename: str, index_info: Optional[Dict[str, Any]]
    ):
        """Set index column information for a file."""
        # In the flat structure, we don't track index info separately
        # since it's embedded in the discrepancies themselves
        pass

    def get_summary(self) -> Dict[str, int]:
        """Generate summary statistics."""
        summary = {
            "total_discrepancies": 0,
            "critical_discrepancies": 0,
            "high_discrepancies": 0,
            "medium_discrepancies": 0,
            "low_discrepancies": 0,
        }

        for discrepancy in self.results:
            summary["total_discrepancies"] += 1
            severity = str(discrepancy["severity"]).lower()
            summary_key = f"{severity}_discrepancies"
            if summary_key in summary:
                summary[summary_key] += 1

        return summary

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary for JSON output."""
        output: Dict[str, Any] = {
            "metadata": {
                "timestamp": self.timestamp,
                "run1_dir": str(self.run1_dir),
                "run2_dir": str(self.run2_dir),
                "versions_run1": self.versions_run1,
                "versions_run2": self.versions_run2,
                "common_versions": self.common_versions,
            },
            "summary": self.get_summary(),
            "results": self.results,
        }
        # Include processing errors if any
        if self.errors:
            output["errors"] = self.errors
        return output

    def to_json(self, indent: int = 2) -> str:
        """Convert report to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)
