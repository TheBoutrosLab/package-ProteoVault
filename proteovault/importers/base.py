"""Source format importer contract and shared parsing operations."""

import math
from abc import ABC, abstractmethod
from pathlib import Path

from proteovault.database import ProteoVaultDatabase


class BaseImporter(ABC):
    """Map one source artifact into the vault without committing transactions."""

    def __init__(
        self,
        path: str | Path,
        database: ProteoVaultDatabase,
        source_id: str,
        stage: str,
        sample_name: str | None = None,
    ):
        self.path = Path(path)
        self.database = database
        self.source_id = source_id
        self.stage = stage
        self.sample_name = sample_name
        self.counts: dict[str, int] = {}

    @abstractmethod
    def import_data(self) -> dict[str, int]:
        """Read, check, and import the source, returning imported record counts."""

    def record_counts(self) -> dict[str, int]:
        for entity, count in self.counts.items():
            self.database.add(
                "import_counts",
                source_id=self.source_id,
                entity=entity,
                row_count=count,
            )
        self.database.flush()
        return self.counts

    @staticmethod
    def number(value) -> float | None:
        if value is None or str(value).strip().lower() in (
            "",
            "na",
            "nan",
            "null",
            "none",
        ):
            return None
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f"Nonfinite numeric value: {value}")
        return result

    @staticmethod
    def run_name(value: str) -> str:
        name = value.replace("\\", "/").rsplit("/", 1)[-1]
        for suffix in (".mzML", ".raw", ".d"):
            if name.lower().endswith(suffix.lower()):
                return name[: -len(suffix)]
        return name
