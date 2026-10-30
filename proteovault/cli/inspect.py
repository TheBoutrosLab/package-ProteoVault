"""Inspect a ProteoVault database."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import logging
from pathlib import Path

from proteovault.database import ProteoVaultDatabase

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

logger = logging.getLogger(__name__)


class VaultInspector:
    """Read database metadata, table sizes, and source provenance."""

    def __init__(self, database: str | Path) -> None:
        self.database = Path(database)

    def summary(self) -> dict:
        logger.info("Inspecting database %s", self.database)
        with ProteoVaultDatabase(self.database, read_only=True) as database:
            summary = {
                "metadata": dict(
                    database.connection.execute(
                        "SELECT key, value FROM metadata"
                    ).fetchall()
                ),
                "analysis": self.records(database, "SELECT * FROM analysis"),
                "table_counts": self.counts(database),
                "sources": self.records(
                    database, "SELECT * FROM source_files ORDER BY relative_path"
                ),
            }

        stage_layers = {
            "quantitation": "core",
            "quant-details": "quant-details",
            "validation": "core",
            "inference": "core",
            "search": "raw-search",
            "MSBooster": "msbooster",
            "qc": "qc",
            "reference": "fasta",
        }
        summary["layers"] = sorted(
            {
                stage_layers[source["stage"]]
                for source in summary["sources"]
                if source["stage"] in stage_layers
            }
        )
        logger.info(
            "Inspection complete: %d tables, %d source files",
            len(summary["table_counts"]),
            len(summary["sources"]),
        )
        return summary

    def format_summary(self, summary: dict) -> str:
        """Render an inspection summary as readable CLI tables."""
        analysis = summary["analysis"]
        details = [
            ("Database", str(self.database)),
            ("Dataset", ", ".join(str(row["dataset_id"]) for row in analysis)),
            ("Pipeline", ", ".join(str(row["pipeline"]) for row in analysis)),
            (
                "Pipeline version",
                ", ".join(str(row["pipeline_version"] or "-") for row in analysis),
            ),
            ("Layers", ", ".join(summary["layers"]) or "-"),
        ]
        details.extend(sorted(summary["metadata"].items()))
        counts = [(name, str(count)) for name, count in summary["table_counts"].items()]
        sources = [
            (
                str(source["stage"]),
                str(source["format"]),
                str(source["size_bytes"]),
                str(source["relative_path"]),
            )
            for source in summary["sources"]
        ]
        return "\n\n".join(
            (
                "Analysis\n" + self.format_table(("Field", "Value"), details),
                "Table counts\n" + self.format_table(("Table", "Rows"), counts),
                "Sources\n"
                + self.format_table(
                    ("Stage", "Format", "Bytes", "Relative path"), sources
                ),
            )
        )

    @staticmethod
    def format_table(headers: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
        """Align columns in a plain-text table."""
        rows = [tuple(str(value) for value in row) for row in rows]
        widths = [
            max(len(header), *(len(row[index]) for row in rows))
            for index, header in enumerate(headers)
        ]
        lines = [
            "  " + "  ".join(header.ljust(width) for header, width in zip(headers, widths)),
            "  " + "  ".join("-" * width for width in widths),
        ]
        lines.extend(
            "  " + "  ".join(value.ljust(width) for value, width in zip(row, widths))
            for row in rows
        )
        return "\n".join(lines)

    def table_counts(self) -> dict[str, int]:
        with ProteoVaultDatabase(self.database, read_only=True) as database:
            return self.counts(database)

    def sources(self) -> list[dict]:
        with ProteoVaultDatabase(self.database, read_only=True) as database:
            return self.records(
                database, "SELECT * FROM source_files ORDER BY relative_path"
            )

    @staticmethod
    def records(database: ProteoVaultDatabase, query: str) -> list[dict]:
        result = database.connection.execute(query)
        names = [column[0] for column in result.description]
        return [dict(zip(names, row)) for row in result.fetchall()]

    @staticmethod
    def counts(database: ProteoVaultDatabase) -> dict[str, int]:
        tables = database.connection.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main' AND table_type='BASE TABLE' ORDER BY table_name"
        ).fetchall()
        return {
            table: database.connection.execute(
                'SELECT count(*) FROM "' + table.replace('"', '""') + '"'
            ).fetchone()[0]
            for (table,) in tables
        }


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "inspect", help="Show metadata, source files, and record counts"
    )
    parser.add_argument("database", type=Path)
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    inspector = VaultInspector(args.database)
    summary = inspector.summary()
    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(inspector.format_summary(summary))
    return 0
