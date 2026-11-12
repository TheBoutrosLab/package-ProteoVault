"""Resolve acquisition filenames to metapipeline sample names."""

import csv
from pathlib import Path


class SampleManifest:
    """Map the published DIA-NN input paths to canonical sample names."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.by_path = {}
        self.by_stem = {}
        self.samples = set()
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            path_columns = set(fields) & {"path", "mzML"}
            if (len(fields) != len(set(fields)) or "sample_name" not in fields
                    or len(path_columns) != 1):
                raise ValueError(
                    f"Unexpected quantitation manifest columns in {self.path}; "
                    "require sample_name and exactly one of path or mzML "
                    "without duplicate headers"
                )
            path_column = path_columns.pop()
            fractions = set()
            for line, row in enumerate(reader, start=2):
                sample = (row["sample_name"] or "").strip()
                acquisition = (row[path_column] or "").strip()
                if not sample or not acquisition:
                    raise ValueError(
                        f"Incomplete quantitation manifest row: {self.path}:{line}"
                    )
                if "fraction_id" in fields:
                    fraction = (row["fraction_id"] or "").strip()
                    if not fraction or (sample, fraction) in fractions:
                        raise ValueError(
                            f"Missing or duplicate fraction for {sample}: "
                            f"{self.path}:{line}"
                        )
                    fractions.add((sample, fraction))
                normalized = acquisition.replace("\\", "/")
                stem = Path(normalized).name
                for suffix in (".mzML", ".raw", ".d"):
                    if stem.lower().endswith(suffix.lower()):
                        stem = stem[: -len(suffix)]
                        break
                if normalized in self.by_path and self.by_path[normalized] != sample:
                    raise ValueError(
                        f"Conflicting sample names for {acquisition}: {self.path}"
                    )
                self.by_path[normalized] = sample
                self.by_stem.setdefault(stem, set()).add(sample)
                self.samples.add(sample)
        if not self.samples:
            raise ValueError(f"Empty quantitation manifest: {self.path}")

    def sample_name(self, acquisition: str) -> str:
        """Resolve a DIA-NN path or run name without guessing ambiguous stems."""
        normalized = acquisition.replace("\\", "/")
        if normalized in self.by_path:
            return self.by_path[normalized]
        if normalized in self.samples:
            return normalized
        stem = Path(normalized).name
        for suffix in (".mzML", ".raw", ".d"):
            if stem.lower().endswith(suffix.lower()):
                stem = stem[: -len(suffix)]
                break
        matches = self.by_stem.get(stem, set())
        if len(matches) == 1:
            return next(iter(matches))
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous sample mapping for acquisition {acquisition!r}"
            )
        if stem in self.samples:
            return stem
        raise ValueError(f"Acquisition {acquisition!r} is absent from {self.path}")
