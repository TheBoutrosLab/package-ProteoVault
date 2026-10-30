"""Import reference protein sequences without rewriting accession prefixes."""

import gzip
import logging

from proteovault.importers.base import BaseImporter

logger = logging.getLogger(__name__)


class FASTAImporter(BaseImporter):
    """Load FASTA entries using their original first-token accessions."""

    def __init__(self, *args, only_existing: bool = False, **kwargs):
        super().__init__(*args, **kwargs)
        self.only_existing = only_existing

    def import_data(self):
        self.counts = {"protein_entries": 0}
        self.seen_accessions = set()
        self.existing = {}
        if self.only_existing:
            self.database.flush()
            self.existing = dict(
                self.database.connection.execute(
                    "SELECT accession, sequence FROM proteins"
                ).fetchall()
            )
        accession = None
        description = None
        sequence = []
        opener = gzip.open if self.path.name.lower().endswith(".gz") else open
        with opener(self.path, "rt", encoding="utf-8-sig") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if accession is not None:
                        self.add_entry(accession, description, sequence)
                    header = line[1:].split(maxsplit=1)
                    if not header:
                        raise ValueError(f"Empty FASTA header: {self.path}")
                    accession = header[0]
                    description = header[1] if len(header) > 1 else None
                    sequence = []
                elif accession is None:
                    raise ValueError(f"Sequence before FASTA header: {self.path}")
                elif not self.only_existing or accession in self.existing:
                    sequence.append(line)
        if accession is not None:
            self.add_entry(accession, description, sequence)
        if not self.counts["protein_entries"]:
            if self.only_existing:
                raise ValueError(
                    f"No FASTA accessions match proteins in the vault: {self.path}. "
                    "Check that this is the search FASTA and that header identifiers match."
                )
            raise ValueError(f"Empty FASTA: {self.path}")
        if self.only_existing:
            logger.info(
                "FASTA provided sequences for %d of %d identified proteins",
                self.counts["protein_entries"],
                len(self.existing),
            )
        return self.record_counts()

    def add_entry(self, accession, description, sequence):
        if self.only_existing and accession not in self.existing:
            return
        if accession in self.seen_accessions:
            raise ValueError(f"Duplicate FASTA accession {accession!r}: {self.path}")
        joined = "".join(sequence)
        if not joined:
            raise ValueError(f"Empty FASTA sequence for {accession}")
        existing_sequence = self.existing.get(accession)
        if existing_sequence and existing_sequence != joined:
            raise ValueError(
                f"FASTA sequence conflicts with existing sequence for {accession!r}"
            )
        self.database.protein(accession, joined, description)
        self.seen_accessions.add(accession)
        self.counts["protein_entries"] += 1
        if self.counts["protein_entries"] % 2000 == 0:
            self.database.flush()
