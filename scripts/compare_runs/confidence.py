
from enum import Enum

class Confidence(Enum):
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"

    def __str__(self) -> str:
        return self.value
