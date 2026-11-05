"""Import the published cohort QC table without changing metric semantics."""

import csv

from proteovault.importers.base import BaseImporter


class QCImporter(BaseImporter):
    """Stream one merged qc-ProteomicDIA metrics TSV into the vault."""

    FIELDS = [
        "dataset_id",
        "layer",
        "source_process",
        "sample_name",
        "search_engine",
        "fraction_id",
        "metric_name",
        "group",
        "value",
        "unit",
        "status",
        "source_file",
    ]

    def import_data(self) -> dict[str, int]:
        dataset = self.database.connection.execute(
            "SELECT dataset_id FROM analysis"
        ).fetchall()[0][0]
        self.database.flush()
        quantified_samples = {
            name for (name,) in self.database.connection.execute(
                "SELECT DISTINCT s.sample_name FROM quant_channels c "
                "JOIN samples s USING(sample_id)"
            ).fetchall()
        }
        self.counts = {"qc_metrics": 0}
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.FIELDS:
                raise ValueError(
                    f"Unexpected cohort QC columns in {self.path}: {reader.fieldnames}"
                )
            for row_number, row in enumerate(reader, start=1):
                if row["dataset_id"] != dataset:
                    raise ValueError(
                        f"{self.path}: row {row_number} belongs to "
                        f"{row['dataset_id']!r}, expected {dataset!r}"
                    )
                if not row["layer"] or not row["metric_name"]:
                    raise ValueError(
                        f"{self.path}: row {row_number} lacks layer or metric_name"
                    )
                if row["status"] not in ("available", "unavailable"):
                    raise ValueError(
                        f"{self.path}: row {row_number} has invalid status "
                        f"{row['status']!r}"
                    )
                value = row["value"] or None
                if row["status"] == "available" and value is None:
                    raise ValueError(
                        f"{self.path}: row {row_number} has no available value"
                    )
                sample_name = (
                    self.database.sample_manifest.sample_name(row["sample_name"])
                    if row["sample_name"] and self.database.sample_manifest
                    else row["sample_name"]
                )
                if sample_name and sample_name not in quantified_samples:
                    raise ValueError(
                        f"{self.path}: row {row_number} QC sample {row['sample_name']!r} "
                        "does not match a quantified sample; provide a matching "
                        "manifests/quantitation.csv"
                    )
                self.database.add(
                    "qc_metrics",
                    source_id=self.source_id,
                    row_number=row_number,
                    dataset_id=row["dataset_id"],
                    layer=row["layer"],
                    source_process=row["source_process"],
                    sample_name=sample_name,
                    search_engine=row["search_engine"],
                    fraction_id=row["fraction_id"],
                    metric_name=row["metric_name"],
                    group=row["group"],
                    value=value,
                    unit=row["unit"],
                    status=row["status"],
                    source_file=row["source_file"],
                )
                self.counts["qc_metrics"] += 1
                if row_number % 5000 == 0:
                    self.database.flush()
        return self.record_counts()
