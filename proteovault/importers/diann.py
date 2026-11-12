"""Import DIA-NN matrices and long reports without implicit confidence filtering."""

import csv

from proteovault.importers.base import BaseImporter


class DIANNImporter(BaseImporter):
    """Keep report and matrix quantities distinct and retain their measurement methods."""

    def import_data(self) -> dict[str, int]:
        self.counts = {"input_rows": 0, "protein_quant": 0, "precursor_quant": 0}
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)):
                raise ValueError(f"Duplicate DIA-NN columns: {self.path}")
            if "Protein.Group" not in fields:
                raise ValueError(f"Missing Protein.Group column: {self.path}")
            is_matrix = self.path.name.endswith((".pg_matrix.tsv", ".pr_matrix.tsv"))
            precursor = self.path.name.endswith(".pr_matrix.tsv") or not is_matrix
            if precursor:
                required = {
                    "Modified.Sequence",
                    "Stripped.Sequence",
                    "Precursor.Charge",
                }
                if not required.issubset(fields):
                    raise ValueError(
                        f"Missing precursor columns in {self.path}: {sorted(required - set(fields))}"
                    )
            channels = [
                f for f in fields
                if f.lower().endswith((".raw", ".mzml", ".d"))
                or (is_matrix and self.database.sample_manifest
                    and f in self.database.sample_manifest.samples)
            ]
            if is_matrix and not channels:
                raise ValueError(f"No sample columns in {self.path}")
            if not is_matrix and "Run" not in fields:
                raise ValueError(f"Missing Run in {self.path}")
            if not is_matrix and precursor and "Precursor.Quantity" not in fields:
                raise ValueError(f"Missing Precursor.Quantity in {self.path}")
            protein_methods = tuple(
                name for name in ("PG.Quantity", "PG.Normalised", "PG.MaxLFQ", "PG.TopN")
                if not is_matrix and self.stage == "quantitation" and name in fields
            )
            if is_matrix:
                table = "precursor_quant" if precursor else "protein_quant"
                self.database.register_quant_method(self.source_id, table, "matrix")
            else:
                if precursor:
                    for name in ("Precursor.Quantity", "Precursor.Normalised",
                                 "Precursor.Translated"):
                        if name in fields:
                            self.database.register_quant_method(
                                self.source_id, "precursor_quant", name)
                for name in protein_methods:
                    self.database.register_quant_method(
                        self.source_id, "protein_quant", name)
            protein_values = {}
            seen = set()
            self.groups = {}
            self.precursors = {}
            self.mappings = set()
            self.channel_ids = (
                {
                    column: self.database.channel(
                        self.run_name(column),
                        column if column.lower().endswith((".raw", ".mzml", ".d"))
                        else None,
                    )
                    for column in channels
                }
                if is_matrix
                else {}
            )
            for line, row in enumerate(reader, start=2):
                self.counts["input_rows"] += 1
                group_text = row["Protein.Group"]
                if not group_text:
                    raise ValueError(f"{self.path}:{line}: empty Protein.Group")
                if not is_matrix and not precursor and not protein_methods:
                    continue
                group = self.group(group_text)
                precursor_id = None
                if precursor:
                    precursor_id = self.precursor(row)
                row_key = precursor_id if precursor else group
                if is_matrix:
                    if row_key in seen:
                        raise ValueError(f"{self.path}:{line}: repeated matrix entity")
                    seen.add(row_key)
                    for column in channels:
                        channel = self.channel_ids[column]
                        value = self.number(row[column])
                        if precursor:
                            self.add_precursor(precursor_id, channel, "matrix", value)
                        else:
                            self.database.add(
                                "protein_quant",
                                group_id=group,
                                channel_id=channel,
                                source_id=self.source_id,
                                method="matrix",
                                quantity=value,
                            )
                            self.counts["protein_quant"] += 1
                    if precursor and any(
                        self.number(row[column]) is not None for column in channels
                    ):
                        self.map_precursor(precursor_id, group)
                else:
                    run = row["Run"]
                    if run not in self.channel_ids:
                        self.channel_ids[run] = self.database.channel(
                            self.run_name(run), row.get("File.Name")
                        )
                    channel = self.channel_ids[run]
                    if precursor:
                        key = (precursor_id, channel)
                        if key in seen:
                            raise ValueError(
                                f"{self.path}:{line}: repeated precursor/run; refusing ambiguous quantities"
                            )
                        seen.add(key)
                        for column in (
                            "Precursor.Quantity",
                            "Precursor.Normalised",
                            "Precursor.Translated",
                        ):
                            if column in fields:
                                self.add_precursor(
                                    precursor_id,
                                    channel,
                                    column,
                                    self.number(row[column]),
                                    self.number(row.get("Q.Value")),
                                    self.number(row.get("PEP")),
                                    rt_seconds=self.minutes_to_seconds(row.get("RT")),
                                    rt_start_seconds=self.minutes_to_seconds(row.get("RT.Start")),
                                    rt_stop_seconds=self.minutes_to_seconds(row.get("RT.Stop")),
                                )
                    for method in protein_methods:
                        key = (group, channel, method)
                        value = self.number(row[method])
                        if key not in protein_values:
                            protein_values[key] = value
                        elif value is not None:
                            previous = protein_values[key]
                            if previous is not None and previous != value:
                                raise ValueError(
                                    f"{self.path}:{line}: conflicting {method} "
                                    f"for protein group {group_text!r} and run {run!r}"
                                )
                            protein_values[key] = value
                    if precursor and any(
                        self.number(row.get(column)) is not None
                        for column in (
                            "Precursor.Quantity",
                            "Precursor.Normalised",
                            "Precursor.Translated",
                        )
                    ):
                        self.map_precursor(precursor_id, group)
                if line % 5000 == 0:
                    self.database.flush()
            for (group, channel, method), value in protein_values.items():
                self.database.add(
                    "protein_quant",
                    group_id=group,
                    channel_id=channel,
                    source_id=self.source_id,
                    method=method,
                    quantity=value,
                )
                self.counts["protein_quant"] += 1
                if self.counts["protein_quant"] % 5000 == 0:
                    self.database.flush()
        return self.record_counts()

    def group(self, text):
        if text not in self.groups:
            members = [accession for accession in text.split(";") if accession]
            group = self.database.group(self.source_id, text, members, "DIA-NN")
            if len(self.groups) >= 100000:
                self.groups.clear()
            self.groups[text] = group
        return self.groups[text]

    def precursor(self, row):
        key = (
            row["Stripped.Sequence"],
            row["Modified.Sequence"],
            int(row["Precursor.Charge"]),
        )
        if key not in self.precursors:
            precursor = self.database.precursor(*key, "DIA-NN")
            if len(self.precursors) >= 100000:
                self.precursors.clear()
            self.precursors[key] = precursor
        return self.precursors[key]

    def map_precursor(self, precursor, group):
        key = (precursor, group)
        if key in self.mappings:
            return
        self.database.add("precursor_groups", precursor_id=precursor, group_id=group)
        if len(self.mappings) >= 100000:
            self.mappings.clear()
        self.mappings.add(key)

    def add_precursor(self, precursor, channel, method, value, q_value=None, pep=None,
                      rt_seconds=None, rt_start_seconds=None, rt_stop_seconds=None):
        self.database.add(
            "precursor_quant",
            quant_id=self.database.identity("precursor_quant", self.source_id,
                                            precursor, channel, method),
            precursor_id=precursor,
            channel_id=channel,
            source_id=self.source_id,
            method=method,
            quantity=value,
            q_value=q_value,
            pep=pep,
            rt_seconds=rt_seconds,
            rt_start_seconds=rt_start_seconds,
            rt_stop_seconds=rt_stop_seconds,
        )
        self.counts["precursor_quant"] += 1

    def minutes_to_seconds(self, value):
        minutes = self.number(value)
        return None if minutes is None else minutes * 60
