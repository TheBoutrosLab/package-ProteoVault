"""Import UniProt accession cross-references for proteins already in a vault."""

import gzip

from proteovault.importers.base import BaseImporter


def uniprot_accession(identifier: str) -> str | None:
    """Extract an accession only from an explicit UniProt FASTA identifier."""
    parts = identifier.split("|")
    if (len(parts) != 3 or parts[0] not in ("sp", "tr") or
            not parts[1] or not parts[2] or any(" " in part for part in parts)):
        return None
    return parts[1]


class UniProtIdMappingImporter(BaseImporter):
    """Stream a three-column UniProt mapping archive into protein identifiers."""

    def import_data(self) -> dict[str, int]:
        self.database.flush()
        proteins = self.database.connection.execute(
            "SELECT protein_id, accession, is_decoy, is_contaminant FROM proteins"
        ).fetchall()
        by_accession = {}
        skipped = 0
        for protein_id, identifier, is_decoy, is_contaminant in proteins:
            accession = None if is_decoy or is_contaminant else uniprot_accession(identifier)
            if accession is None:
                skipped += 1
                continue
            if accession in by_accession and by_accession[accession] != protein_id:
                raise ValueError(f"Ambiguous UniProt accession in vault: {accession}")
            by_accession[accession] = protein_id
        if not by_accession:
            raise ValueError("Vault has no unambiguous UniProt protein accessions")

        self.counts = {
            "input_rows": 0,
            "protein_identifiers": 0,
            "matched_proteins": 0,
            "unmatched_proteins": 0,
            "skipped_proteins": skipped,
        }
        matched = set()
        seen = set()
        with gzip.open(self.path, "rt", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                self.counts["input_rows"] += 1
                columns = line.rstrip("\r\n").split("\t")
                if len(columns) != 3:
                    raise ValueError(
                        f"{self.path}:{line_number}: expected three tab-separated columns"
                    )
                accession, mapping_type, mapped_id = columns
                if not accession or not mapping_type:
                    raise ValueError(f"{self.path}:{line_number}: missing accession or mapping type")
                protein_id = by_accession.get(accession)
                if protein_id is None or not mapped_id or mapped_id == "-":
                    continue
                matched.add(accession)
                key = (protein_id, mapping_type, mapped_id)
                if key in seen:
                    continue
                seen.add(key)
                self.database.add(
                    "protein_identifiers",
                    protein_id=protein_id,
                    source_id=self.source_id,
                    source_accession=accession,
                    mapping_type=mapping_type,
                    mapped_id=mapped_id,
                )
                self.counts["protein_identifiers"] += 1
                if self.counts["protein_identifiers"] % 5000 == 0:
                    self.database.flush()
        if not matched:
            raise ValueError(f"No vault UniProt accessions matched {self.path}")
        self.counts["matched_proteins"] = len(matched)
        self.counts["unmatched_proteins"] = len(by_accession) - len(matched)
        return self.record_counts()
