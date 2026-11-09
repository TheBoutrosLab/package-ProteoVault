"""Validate the schema, relationships, and imported record counts."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import logging
from pathlib import Path

import duckdb

from proteovault.database import ProteoVaultDatabase
from proteovault.schema import SCHEMA_VERSION

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

logger = logging.getLogger(__name__)


class VaultValidator:
    """Return machine-readable validation findings without modifying the database."""

    def __init__(self, database: str | Path) -> None:
        self.database = Path(database)

    def validate(self) -> dict:
        logger.info("Validating database %s", self.database)
        errors = []
        warnings = []
        try:
            with ProteoVaultDatabase(self.database, read_only=True) as database:
                connection = database.connection
                required = self.expected_columns()
                actual = {}
                for table, column in connection.execute(
                    "SELECT table_name, column_name FROM information_schema.columns WHERE table_schema='main'"
                ).fetchall():
                    actual.setdefault(table, set()).add(column)
                for table, columns in required.items():
                    missing = columns - actual.get(table, set())
                    if missing:
                        errors.append(f"{table}: missing columns {sorted(missing)}")
                if errors:
                    return self.report(errors, warnings)
                metadata = dict(
                    connection.execute("SELECT key,value FROM metadata").fetchall()
                )
                if metadata.get("schema_version") != SCHEMA_VERSION:
                    errors.append(
                        f"Unsupported schema version: {metadata.get('schema_version')}"
                    )
                if metadata.get("format") != "ProteoVault":
                    errors.append("Missing ProteoVault format identity")
                if (
                    connection.execute("SELECT count(*) FROM analysis").fetchone()[0]
                    != 1
                ):
                    errors.append("Expected exactly one analysis")
                for table, key, parent, parent_key in self.relationships():
                    query = (
                        f'SELECT count(*) FROM "{table}" t LEFT JOIN "{parent}" p '
                        f'ON t."{key}"=p."{parent_key}" WHERE t."{key}" IS NOT NULL '
                        f'AND p."{parent_key}" IS NULL'
                    )
                    count = connection.execute(query).fetchone()[0]
                    if count:
                        errors.append(f"{table}.{key}: {count} orphan references")
                for table, key in [
                    ("samples", "sample_id"),
                    ("runs", "run_id"),
                    ("psms", "psm_id"),
                    ("protein_groups", "group_id"),
                    ("spectra", "spectrum_id"),
                    ("peptides", "peptide_id"),
                    ("peptidoforms", "peptidoform_id"),
                    ("precursors", "precursor_id"),
                    ("precursor_quant", "quant_id"),
                    ("proteins", "protein_id"),
                    ("quant_channels", "channel_id"),
                    ("tmt_plexes", "plex_id"),
                    ("tmt_channels", "channel_id"),
                    ("tmt_runs", "run_id"),
                    ("peptide_quant_features", "feature_id"),
                ]:
                    total, unique = connection.execute(
                        f'SELECT count(*),count(DISTINCT "{key}") FROM "{table}"'
                    ).fetchone()
                    if total != unique:
                        errors.append(f"Duplicate or null keys in {table}")
                for table, columns in (
                    ("source_files", "path, format, stage"),
                    ("samples", "sample_name"),
                    ("tmt_plexes", "plex_name"),
                    ("tmt_channels", "plex_id, reporter_label"),
                    ("runs", "run_name"),
                    ("spectra", "run_id, native_id"),
                    ("peptides", "sequence"),
                    ("peptidoforms", "modified_sequence, notation"),
                    ("precursors", "peptidoform_id, charge"),
                    ("proteins", "accession"),
                    ("quant_methods", "name"),
                    ("grouping_methods", "name"),
                    ("processing_stages", "name"),
                    ("score_types", "name"),
                ):
                    duplicate = connection.execute(
                        f"SELECT 1 FROM {table} GROUP BY {columns} HAVING count(*) > 1 LIMIT 1"
                    ).fetchone()
                    if duplicate:
                        errors.append(f"{table}: duplicate natural key ({columns})")
                wrong_reference = connection.execute(
                    "SELECT count(*) FROM tmt_plexes p "
                    "LEFT JOIN tmt_channels c ON p.reference_channel_id=c.channel_id "
                    "WHERE p.reference_channel_id IS NOT NULL "
                    "AND (c.channel_id IS NULL OR c.plex_id<>p.plex_id)"
                ).fetchall()[0][0]
                if wrong_reference:
                    errors.append(
                        f"tmt_plexes: {wrong_reference} missing or cross-plex reference channels"
                    )
                mislabeled_channels = connection.execute(
                    "SELECT count(*) FROM tmt_channels t "
                    "JOIN tmt_plexes p USING(plex_id) "
                    "JOIN quant_channels c USING(channel_id) "
                    "JOIN samples s USING(sample_id) "
                    "WHERE s.sample_name <> p.plex_name || '|' || t.reporter_label"
                ).fetchall()[0][0]
                if mislabeled_channels:
                    errors.append(
                        f"tmt_channels: {mislabeled_channels} sample names do not match plex|reporter"
                    )
                invalid_tmt_runs = connection.execute(
                    "SELECT count(*) FROM tmt_runs t JOIN runs r USING(run_id) "
                    "WHERE r.sample_id IS NOT NULL"
                ).fetchall()[0][0]
                if invalid_tmt_runs:
                    errors.append(
                        f"tmt_runs: {invalid_tmt_runs} runs also point to individual samples"
                    )
                for table in ("protein_quant", "precursor_quant", "peptide_quant"):
                    invalid = connection.execute(
                        f"SELECT count(*) FROM {table} WHERE quantity < 0 OR NOT isfinite(quantity)"
                    ).fetchone()[0]
                    if invalid:
                        errors.append(
                            f"{table}: {invalid} negative/nonfinite quantities"
                        )
                invalid = connection.execute(
                    "SELECT count(*) FROM precursor_quant WHERE q_value NOT BETWEEN 0 AND 1 OR pep NOT BETWEEN 0 AND 1"
                ).fetchone()[0]
                if invalid:
                    errors.append(
                        f"precursor_quant: {invalid} invalid confidence values"
                    )
                invalid = connection.execute(
                    "SELECT count(*) FROM precursor_quant q JOIN source_files f USING(source_id) "
                    "WHERE f.format IN ('tmt_consensus_scan','tmt_consensus_feature') "
                    "AND (q.source_element_id IS NULL OR q.source_map_id IS NULL "
                    "OR q.rt_seconds IS NULL "
                    "OR (f.format='tmt_consensus_scan' AND q.spectrum_id IS NULL) "
                    "OR (f.format='tmt_consensus_feature' AND q.spectrum_id IS NOT NULL))"
                ).fetchone()[0]
                if invalid:
                    errors.append(f"precursor_quant: {invalid} invalid TMT observation links")
                for source_id, entity, expected in connection.execute(
                    "SELECT source_id,entity,row_count FROM import_counts "
                    "WHERE entity IN ('psms','protein_quant','precursor_quant','peptide_quant','peptide_quant_features','pin_features','qc_metrics','protein_identifiers')"
                ).fetchall():
                    actual_count = connection.execute(
                        f"SELECT count(*) FROM {entity} WHERE source_id=?", [source_id]
                    ).fetchone()[0]
                    if expected != actual_count:
                        errors.append(
                            f"{entity} for {source_id}: expected {expected}, found {actual_count}"
                        )
                invalid = connection.execute(
                    "SELECT count(*) FROM qc_metrics WHERE "
                    "status='available' AND value IS NULL"
                ).fetchone()[0]
                if invalid:
                    errors.append(f"qc_metrics: {invalid} missing available values")
                if (
                    connection.execute("SELECT count(*) FROM protein_quant").fetchone()[
                        0
                    ]
                    == 0
                ):
                    errors.append("No protein quantities imported")
                if connection.execute("SELECT count(*) FROM psms").fetchone()[0] == 0:
                    warnings.append("No sample-level PSMs imported")
                warnings.extend(json.loads(metadata.get("warnings", "[]")))
        except (duckdb.Error, OSError, ValueError) as error:
            errors.append(str(error))
        return self.report(errors, warnings)

    def report(self, errors: list[str], warnings: list[str]) -> dict:
        warnings = list(dict.fromkeys(warnings))
        for warning in warnings:
            logger.warning("%s", warning)
        for error in errors:
            logger.error("%s", error)
        if errors:
            logger.error(
                "Validation failed for %s: %d errors", self.database, len(errors)
            )
        else:
            logger.info(
                "Validation passed for %s (%d warnings)", self.database, len(warnings)
            )
        return dict(valid=not errors, errors=errors, warnings=warnings)

    @staticmethod
    def expected_columns() -> dict[str, set[str]]:
        from importlib.resources import files

        with duckdb.connect() as connection:
            connection.execute(
                files("proteovault.schema").joinpath("v1.sql").read_text()
            )
            result = {}
            for table, column in connection.execute(
                "SELECT table_name,column_name FROM information_schema.columns WHERE table_schema='main'"
            ).fetchall():
                result.setdefault(table, set()).add(column)
            return result

    @staticmethod
    def relationships() -> list[tuple[str, str, str, str]]:
        return [
            ("runs", "sample_id", "samples", "sample_id"),
            ("spectra", "run_id", "runs", "run_id"),
            ("spectra", "source_id", "source_files", "source_id"),
            ("peptidoforms", "peptide_id", "peptides", "peptide_id"),
            ("precursors", "peptidoform_id", "peptidoforms", "peptidoform_id"),
            ("psms", "spectrum_id", "spectra", "spectrum_id"),
            ("psms", "precursor_id", "precursors", "precursor_id"),
            ("psms", "source_id", "source_files", "source_id"),
            ("psm_scores", "psm_id", "psms", "psm_id"),
            ("psm_scores", "score_id", "score_types", "score_id"),
            ("psms", "stage_id", "processing_stages", "stage_id"),
            ("protein_quant", "method_id", "quant_methods", "method_id"),
            ("protein_groups", "method_id", "grouping_methods", "method_id"),
            ("precursor_quant", "method_id", "quant_methods", "method_id"),
            ("psm_features", "psm_id", "psms", "psm_id"),
            ("psm_metadata", "psm_id", "psms", "psm_id"),
            ("peptide_proteins", "peptide_id", "peptides", "peptide_id"),
            ("peptide_proteins", "protein_id", "proteins", "protein_id"),
            ("psm_proteins", "psm_id", "psms", "psm_id"),
            ("psm_proteins", "protein_id", "proteins", "protein_id"),
            ("protein_groups", "source_id", "source_files", "source_id"),
            ("protein_group_members", "group_id", "protein_groups", "group_id"),
            ("protein_group_members", "protein_id", "proteins", "protein_id"),
            ("protein_identifiers", "protein_id", "proteins", "protein_id"),
            ("protein_identifiers", "source_id", "source_files", "source_id"),
            ("quant_channels", "sample_id", "samples", "sample_id"),
            ("quant_channels", "run_id", "runs", "run_id"),
            ("tmt_plexes", "reference_channel_id", "tmt_channels", "channel_id"),
            ("tmt_channels", "channel_id", "quant_channels", "channel_id"),
            ("tmt_channels", "plex_id", "tmt_plexes", "plex_id"),
            ("tmt_runs", "run_id", "runs", "run_id"),
            ("tmt_runs", "plex_id", "tmt_plexes", "plex_id"),
            ("precursor_quant", "precursor_id", "precursors", "precursor_id"),
            ("precursor_quant", "channel_id", "quant_channels", "channel_id"),
            ("precursor_quant", "source_id", "source_files", "source_id"),
            ("precursor_quant", "spectrum_id", "spectra", "spectrum_id"),
            ("protein_quant", "group_id", "protein_groups", "group_id"),
            ("protein_quant", "channel_id", "quant_channels", "channel_id"),
            ("protein_quant", "source_id", "source_files", "source_id"),
            ("peptide_quant_features", "source_id", "source_files", "source_id"),
            ("peptide_quant_features", "peptide_id", "peptides", "peptide_id"),
            ("peptide_quant", "feature_id", "peptide_quant_features", "feature_id"),
            ("peptide_quant", "channel_id", "quant_channels", "channel_id"),
            ("peptide_quant", "source_id", "source_files", "source_id"),
            ("peptide_quant", "method_id", "quant_methods", "method_id"),
            ("precursor_groups", "precursor_id", "precursors", "precursor_id"),
            ("precursor_groups", "group_id", "protein_groups", "group_id"),
            ("pin_features", "source_id", "source_files", "source_id"),
            ("import_counts", "source_id", "source_files", "source_id"),
            ("qc_metrics", "source_id", "source_files", "source_id"),
        ]


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "validate", help="Check schema, relationships, and counts"
    )
    parser.add_argument("database", type=Path)
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    report = VaultValidator(args.database).validate()
    print(json.dumps(report, indent=2))
    return 0 if report["valid"] else 1
