import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, TypeAlias

from .discrepancy import Discrepancy as ImportedDiscrepancy
from .errors import ComparisonError


Discrepancy: TypeAlias = ImportedDiscrepancy  # type: ignore[assignment]


class ComparisonReport:
    """Collects and manages all discrepancies found during comparison."""

    def __init__(self, run1_dir: Path, run2_dir: Path):
        self.run1_dir = run1_dir
        self.run2_dir = run2_dir
        self.timestamp = datetime.now().isoformat()
        self.versions_run1: List[str] = []
        self.versions_run2: List[str] = []
        self.common_versions: List[str] = []
        self.results: List[Discrepancy] = []
        # Processing errors during comparison (exceptions, etc.)
        self.errors: List[ComparisonError] = []

    def add_discrepancy(self, discrepancy: Discrepancy):
        """Add a discrepancy to the report."""
        self.results.append(discrepancy)

    def add_error(self, error: ComparisonError):
        """Add a processing error to the report."""
        self.errors.append(error)

    def mark_file_identical(self):
        """Mark a file as identical between runs."""
        # In the flat structure, we don't track identical files separately
        # since they don't generate discrepancies
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
            severity = str(discrepancy.severity).lower()
            summary_key = f"{severity}_discrepancies"
            if summary_key in summary:
                summary[summary_key] += 1

        return summary

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to dictionary for JSON output."""
        results = [discrepancy.to_dict() for discrepancy in self.results]
        errors = [error.to_dict() for error in self.errors]
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
            "results": results,
            "errors": errors,
        }
        return output

    def to_json(self, indent: int = 2) -> str:
        """Convert report to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)
