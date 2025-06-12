from dataclasses import dataclass
from typing import Any, Dict, Optional
from typing import Union

from .severity import Severity
from .confidence import Confidence


@dataclass(frozen=True)
class FileReadErrorDiscrepancy:
    severity: Severity
    confidence: Confidence
    description: str
    file: str
    version: str
    run: str
    file_path: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.__class__.__name__,
            "description": self.description,
            "severity": str(self.severity),
            "confidence": str(self.confidence),
            "location": {
                "file": self.file,
                "version": self.version,
                "run": self.run,
                "file_path": self.file_path,
            },
        }


@dataclass(frozen=True)
class EmptyFileDiscrepancy:
    severity: Severity
    confidence: Confidence
    description: str
    file: str
    version: str
    run: str
    file_path: str
    file_size: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.__class__.__name__,
            "description": self.description,
            "severity": str(self.severity),
            "confidence": str(self.confidence),
            "location": {
                "file": self.file,
                "version": self.version,
                "run": self.run,
                "file_path": self.file_path,
                "file_size": self.file_size,
            },
        }


@dataclass(frozen=True)
class NoIndexColumnDiscrepancy:
    severity: Severity
    confidence: Confidence
    description: str
    file: str
    version: str
    run: str
    reason: str
    index_analysis: Optional[Dict[str, Any]]
    has_duplicate_columns_run1: bool
    has_duplicate_columns_run2: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.__class__.__name__,
            "description": self.description,
            "severity": str(self.severity),
            "confidence": str(self.confidence),
            "reason": self.reason,
            "index_analysis": self.index_analysis,
            "has_duplicate_columns_run1": self.has_duplicate_columns_run1,
            "has_duplicate_columns_run2": self.has_duplicate_columns_run2,
            "location": {
                "file": self.file,
                "version": self.version,
                "run": self.run,
            },
        }


# Union type for processing error classes
ComparisonError = Union[
    FileReadErrorDiscrepancy,
    EmptyFileDiscrepancy,
    NoIndexColumnDiscrepancy,
]
