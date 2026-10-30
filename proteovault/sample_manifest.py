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
            if reader.fieldnames != ["sample_name", "mzML"]:
                raise ValueError(
                    f"Unexpected quantitation manifest columns: {self.path}"
                )
            for line, row in enumerate(reader, start=2):
                sample = row["sample_name"].strip()
                mzml = row["mzML"].strip()
                if not sample or not mzml:
                    raise ValueError(
                        f"Incomplete quantitation manifest row: {self.path}:{line}"
                    )
                normalized = mzml.replace("\\", "/")
                stem = Path(normalized).name
                for suffix in (".mzML", ".raw", ".d"):
                    if stem.lower().endswith(suffix.lower()):
                        stem = stem[: -len(suffix)]
                        break
                if normalized in self.by_path and self.by_path[normalized] != sample:
                    raise ValueError(
                        f"Conflicting sample names for {mzml}: {self.path}"
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
