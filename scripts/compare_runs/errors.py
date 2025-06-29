from dataclasses import dataclass
from typing import Any, Dict, List, Optional
from typing import Union

from .severity import Severity
from .confidence import Confidence


@dataclass(frozen=True)
class FileReadError:
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
class EmptyFile:
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
class NoIndexColumn:
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


@dataclass(frozen=True)
class MultipleIndexColumns:
    severity: Severity
    confidence: Confidence
    description: str
    file: str
    version: str
    run: str
    pattern: str
    matching_columns: List[str]
    matching_pairs: List[tuple]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "type": self.__class__.__name__,
            "description": self.description,
            "severity": str(self.severity),
            "confidence": str(self.confidence),
            "pattern": self.pattern,
            "matching_columns": self.matching_columns,
            "matching_pairs": self.matching_pairs,
            "location": {
                "file": self.file,
                "version": self.version,
                "run": self.run,
            },
        }


# Union type for processing error classes
ComparisonError = Union[
    FileReadError,
    EmptyFile,
    NoIndexColumn,
    MultipleIndexColumns,
]
