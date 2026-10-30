"""Preserve PIN features without guessing ambiguous scan-to-PSM links."""

import csv

from proteovault.importers.base import BaseImporter


class PINImporter(BaseImporter):
    """Retain SpecId/ScanNr and all numeric features for later explicit PSM joins."""

    def import_data(self) -> dict[str, int]:
        self.counts = {"pin_rows": 0, "pin_features": 0}
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            fields = next(reader, [])
            if not {"SpecId", "Label", "ScanNr"}.issubset(fields):
                raise ValueError(f"PIN requires SpecId, Label, ScanNr: {self.path}")
            for index, values in enumerate(reader, start=1):
                if not values:
                    continue
                if len(values) < len(fields):
                    raise ValueError(f"Truncated PIN row {index}: {self.path}")
                row = dict(zip(fields, values))
                if row.get("SpecId") == "DefaultDirection":
                    continue
                self.counts["pin_rows"] += 1
                for name, value in row.items():
                    if name in ("SpecId", "Label", "ScanNr", "Peptide", "Proteins"):
                        continue
                    numeric = self.number(value)
                    if numeric is not None:
                        self.database.add(
                            "pin_features",
                            source_id=self.source_id,
                            row_number=index,
                            spectrum_reference=row["SpecId"],
                            scan_number=row["ScanNr"],
                            label=int(row["Label"]),
                            feature_name=name,
                            value=numeric,
                        )
                        self.counts["pin_features"] += 1
                if index % 5000 == 0:
                    self.database.flush()
        return self.record_counts()
