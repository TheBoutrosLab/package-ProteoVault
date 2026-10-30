"""Database access, batch ingestion, and shared entity identities."""

import hashlib
import json
from contextlib import suppress
from importlib.resources import files
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa


class ProteoVaultDatabase:
    """Own a database connection and insert related entities in dependency order."""

    def __init__(self, path: str | Path, read_only: bool = False):
        self.path = Path(path)
        self.connection = duckdb.connect(str(self.path), read_only=read_only)
        self.buffers: dict[str, dict[tuple, dict]] = {}
        self.columns: dict[str, list[str]] = {}
        self.keys: dict[str, list[str]] = {}
        self.column_sets: dict[str, set[str]] = {}
        self.arrow_schemas: dict[str, pa.Schema] = {}
        self.sample_manifest = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self.connection.close()

    def initialize(self):
        self.connection.execute(
            files("proteovault.schema").joinpath("v1.sql").read_text()
        )
        self.load_layout()

    def load_layout(self):
        tables = self.connection.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema = 'main' AND table_type = 'BASE TABLE' ORDER BY table_name"
        ).fetchall()
        for (table,) in tables:
            rows = self.connection.execute(f'DESCRIBE "{table}"').fetchall()
            self.columns[table] = [row[0] for row in rows]
            self.keys[table] = [row[0] for row in rows if row[3] == "PRI"]
            self.column_sets[table] = set(self.columns[table])
            types = {
                "VARCHAR": pa.string(),
                "DOUBLE": pa.float64(),
                "INTEGER": pa.int32(),
                "BIGINT": pa.int64(),
                "BOOLEAN": pa.bool_(),
            }
            self.arrow_schemas[table] = pa.schema(
                [pa.field(row[0], types[row[1]]) for row in rows]
            )
            self.buffers[table] = {}

    @staticmethod
    def identity(kind: str, *parts: Any) -> int:
        """Derive a stable 63-bit relationship key from explicit identity fields."""
        payload = json.dumps([kind, *parts], ensure_ascii=False, separators=(",", ":"))
        return int.from_bytes(
            hashlib.blake2b(payload.encode(), digest_size=8).digest(), "big"
        ) & ((1 << 63) - 1)

    def add(self, table: str, **values):
        if table not in self.columns:
            raise ValueError(f"Unknown table: {table}")
        lookup = {
            "method": ("quant_methods", "method_id"),
            "stage": ("processing_stages", "stage_id"),
            "name": ("score_types", "score_id"),
        }
        if table == "protein_groups" and "method" in values:
            label = values.pop("method")
            identifier = self.identity("grouping_methods", label) & ((1 << 31) - 1)
            self.add("grouping_methods", method_id=identifier, name=label)
            values["method_id"] = identifier
        field = (
            "method"
            if table in ("protein_quant", "precursor_quant")
            else (
                "stage"
                if table == "psms"
                else "name" if table == "psm_scores" else None
            )
        )
        if field is not None and field in values:
            lookup_table, id_field = lookup[field]
            label = values.pop(field)
            identifier = self.identity(lookup_table, label) & ((1 << 31) - 1)
            self.add(lookup_table, **{id_field: identifier, "name": label})
            values[id_field] = identifier
        unknown = values.keys() - self.column_sets[table]
        if unknown:
            raise ValueError(f"Unknown columns for {table}: {sorted(unknown)}")
        key = tuple(values[name] for name in self.keys[table])
        existing = self.buffers[table].get(key)
        if existing and table in ("proteins", "spectra", "runs"):
            values = {
                **existing,
                **{name: value for name, value in values.items() if value is not None},
            }
        self.buffers[table][key] = values

    def flush(self):
        """Insert bounded Arrow batches in schema dependency order."""
        order = [
            "metadata",
            "analysis",
            "source_files",
            "samples",
            "runs",
            "spectra",
            "peptides",
            "peptidoforms",
            "precursors",
            "processing_stages",
            "score_types",
            "score_aliases",
            "quant_methods",
            "grouping_methods",
            "psms",
            "psm_scores",
            "psm_features",
            "psm_metadata",
            "proteins",
            "peptide_proteins",
            "psm_proteins",
            "protein_groups",
            "protein_group_members",
            "quant_channels",
            "precursor_quant",
            "protein_quant",
            "precursor_groups",
            "pin_features",
            "qc_metrics",
            "import_counts",
        ]
        for table in order:
            rows = self.buffers[table]
            if not rows:
                continue
            arrow = pa.Table.from_pylist(
                list(rows.values()), schema=self.arrow_schemas[table]
            )
            self.connection.register("incoming_rows", arrow)
            try:
                if table == "proteins":
                    self.connection.execute(
                        "INSERT INTO proteins BY NAME SELECT * FROM incoming_rows ON CONFLICT (protein_id) DO UPDATE SET "
                        "sequence=coalesce(excluded.sequence,proteins.sequence), "
                        "description=coalesce(excluded.description,proteins.description)"
                    )
                else:
                    self.connection.execute(
                        f'INSERT INTO "{table}" BY NAME SELECT * FROM incoming_rows ON CONFLICT DO NOTHING'
                    )
            finally:
                with suppress(duckdb.Error):
                    self.connection.unregister("incoming_rows")
            rows.clear()

    def sample(self, name: str) -> str:
        key = self.identity("sample", name)
        self.add("samples", sample_id=key, sample_name=name)
        return key

    def run(self, name: str, file_name: str | None = None) -> str:
        key = self.identity("run", name)
        self.add(
            "runs",
            run_id=key,
            run_name=name,
            file_name=file_name,
            sample_id=self.sample(
                self.sample_manifest.sample_name(file_name or name)
                if self.sample_manifest
                else name
            ),
        )
        return key

    def resolve_run(self, sample_name: str) -> str:
        """Use an exact or unique bounded sample match; never guess between runs."""
        if self.sample_manifest:
            sample_name = self.sample_manifest.sample_name(sample_name)
        import re

        self.flush()
        rows = self.connection.execute(
            "SELECT run_id, run_name, sample_id FROM runs"
        ).fetchall()
        if self.sample_manifest:
            linked = [
                row for row in rows if row[2] == self.identity("sample", sample_name)
            ]
            if len(linked) == 1:
                return linked[0][0]
            if len(linked) > 1:
                raise ValueError(f"Multiple acquisition runs for {sample_name}")
        exact = [row for row in rows if row[1] == sample_name]
        if exact:
            return exact[0][0]
        matches = [
            row
            for row in rows
            if re.search(r"(?:^|_)" + re.escape(sample_name) + r"(?:_|$)", row[1])
        ]
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous acquisition runs for {sample_name}; explicit mapping is required"
            )
        if matches:
            return matches[0][0]
        return self.run(sample_name)

    def channel(self, name: str, file_name: str | None = None) -> str:
        key = self.identity("channel", name)
        self.add(
            "quant_channels",
            channel_id=key,
            sample_id=self.sample(
                self.sample_manifest.sample_name(file_name or name)
                if self.sample_manifest
                else name
            ),
            run_id=self.run(name, file_name),
            label="LFQ",
        )
        return key

    def protein(
        self,
        accession: str,
        sequence: str | None = None,
        description: str | None = None,
    ) -> str:
        key = self.identity("protein", accession)
        self.add(
            "proteins",
            protein_id=key,
            accession=accession,
            sequence=sequence,
            description=description,
            is_decoy=accession.startswith(("DECOY_", "rev_")),
            is_contaminant="contam" in accession.lower(),
        )
        return key

    def precursor(
        self, sequence: str, modified: str, charge: int, notation: str
    ) -> str:
        if not sequence or charge <= 0:
            raise ValueError(f"Invalid precursor: {modified!r}, charge {charge}")
        peptide = self.identity("peptide", sequence)
        form = self.identity("form", modified, notation)
        key = self.identity("precursor", form, charge)
        self.add("peptides", peptide_id=peptide, sequence=sequence)
        self.add(
            "peptidoforms",
            peptidoform_id=form,
            peptide_id=peptide,
            modified_sequence=modified,
            notation=notation,
        )
        self.add("precursors", precursor_id=key, peptidoform_id=form, charge=charge)
        return key

    def group(
        self,
        source_id: str,
        identifier: str,
        accessions: list[str],
        method: str,
        group_type: str = "quantification",
        score: float | None = None,
        score_type: str | None = None,
    ) -> str:
        key = self.identity("group", source_id, identifier, group_type)
        self.add(
            "protein_groups",
            group_id=key,
            source_id=source_id,
            source_group_id=identifier,
            method=method,
            group_type=group_type,
            score=score,
            score_type=score_type,
        )
        for accession in accessions:
            self.add(
                "protein_group_members",
                group_id=key,
                protein_id=self.protein(accession),
            )
        return key
