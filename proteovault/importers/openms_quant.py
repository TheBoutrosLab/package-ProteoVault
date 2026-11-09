"""Import OpenMS ProteinQuantifier protein and peptide abundance tables."""

import csv
import re

from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class OpenMSQuantImporter(BaseImporter):
    """Preserve source rows, abundance columns, and missing values."""

    COMMENT_COLUMN = re.compile(r"(\d+): '([^']+)'(?: \('([^']+)'\))?")

    def import_data(self) -> dict[str, int]:
        protein_level = self.path.name.endswith("_ProteinQuantifier_protein.csv")
        entity = "protein_quant" if protein_level else "peptide_quant"
        self.counts = {entity: 0}
        if not protein_level:
            self.counts["peptide_quant_features"] = 0
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            comments = []

            def records():
                for line in handle:
                    if line.startswith("#"):
                        comments.append(line)
                    elif line.strip():
                        yield line

            reader = csv.DictReader(records(), delimiter="\t")
            fields = reader.fieldnames or []
            feature_name = "protein" if protein_level else "peptide"
            if feature_name not in fields:
                raise ValueError(f"Missing {feature_name} column: {self.path}")
            abundance = [name for name in fields if name.startswith((
                "abundance_", "abundance|"
            ))]
            if not abundance:
                raise ValueError(f"No ProteinQuantifier abundance columns: {self.path}")
            column_map = self.column_map(comments, abundance)
            for row_number, row in enumerate(reader, start=1):
                source_value = row[feature_name]
                if not source_value:
                    raise ValueError(f"Empty {feature_name}: {self.path}:{row_number}")
                if protein_level:
                    members = source_value.split("/")
                    score = self.number(row.get("protein_score"))
                    group = self.database.group(
                        self.source_id, source_value, members,
                        "ProteinQuantifier", "quantification", score,
                        "protein_score" if score is not None else None,
                        int(row["n_proteins"]) if row.get("n_proteins") else None,
                        int(row["n_peptides"]) if row.get("n_peptides") else None,
                    )
                else:
                    sequence = IdXMLImporter.unmodified(source_value)
                    if not sequence:
                        raise ValueError(
                            f"Empty unmodified peptide: {self.path}:{row_number}"
                        )
                    peptide = self.database.identity("peptide", sequence)
                    self.database.add("peptides", peptide_id=peptide,
                                      sequence=sequence)
                    feature = self.database.identity(
                        "openms_peptide_quant", self.source_id, row_number
                    )
                    self.database.add(
                        "peptide_quant_features", feature_id=feature,
                        source_id=self.source_id, source_row_number=row_number,
                        peptide_id=peptide, source_peptide=source_value,
                        source_protein=row.get("protein"),
                        reported_n_proteins=(int(row["n_proteins"])
                                             if row.get("n_proteins") else None),
                        charge=int(row["charge"]) if row.get("charge") else None,
                        fraction=row.get("fraction"),
                    )
                    self.counts["peptide_quant_features"] += 1
                for column, channel in column_map.items():
                    values = dict(
                        channel_id=channel, source_id=self.source_id,
                        method="ProteinQuantifier", quantity=self.number(row[column]),
                    )
                    if protein_level:
                        self.database.add("protein_quant", group_id=group, **values)
                    else:
                        self.database.add("peptide_quant", feature_id=feature, **values)
                    self.counts[entity] += 1
                if row_number % 5000 == 0:
                    self.database.flush()
        return self.record_counts()

    def column_map(self, comments: list[str], abundance: list[str]) -> dict[str, int]:
        """Resolve each abundance column to one registered LFQ or TMT channel."""
        annotations = {}
        for comment in comments:
            if "Files/samples associated with abundance values below:" in comment:
                annotations.update({int(number): (sample, label) for number, sample, label
                                    in self.COMMENT_COLUMN.findall(comment)})
        mapped = {}
        for column in abundance:
            number = re.fullmatch(r"abundance_(\d+)", column)
            detailed = re.fullmatch(r"abundance\|(.+)\|ch(\d+)", column)
            numbered_sample = re.fullmatch(
                r"abundance_sample(\d+)\[[^]]*\]", column
            )
            if number:
                entry = annotations.get(int(number.group(1)))
                if entry is None:
                    raise ValueError(f"Missing abundance column annotation: {column}")
                sample, label = entry
                if self.database.dda_mode == "TMT":
                    reporter = label.rsplit("_", 1)[-1]
                    key = (sample, reporter)
                else:
                    key = sample
            elif detailed:
                sample, index = detailed.group(1), int(detailed.group(2))
                if self.database.dda_mode == "TMT":
                    labels = self.database.dda_reporters.get(sample)
                    if not labels or index < 1 or index > len(labels):
                        raise ValueError(f"Unmapped TMT column: {column}")
                    key = (sample, labels[index - 1])
                else:
                    key = sample
            elif numbered_sample and self.database.dda_mode == "TMT":
                if len(self.database.dda_reporters) != 1:
                    raise ValueError(
                        f"Cannot map {column} across multiple TMT plexes: "
                        "numbered protein columns do not identify a plex or reporter; "
                        "export protein quantities with explicit file and channel headers"
                    )
                plex, labels = next(iter(self.database.dda_reporters.items()))
                index = int(numbered_sample.group(1))
                if index < 1 or index > len(labels):
                    raise ValueError(f"Unmapped TMT column: {column}")
                key = (plex, labels[index - 1])
            else:
                raise ValueError(f"Unsupported abundance column: {column}")
            try:
                mapped[column] = self.database.dda_channels[key]
            except KeyError as error:
                raise ValueError(f"No DDA channel mapping for {column}: {key}") from error
        return mapped
