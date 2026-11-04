"""Repair source provenance after moving a ProteoVault analysis."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from argparse import ArgumentParser, Namespace, _SubParsersAction
from pathlib import Path

from proteovault.cli.validate import VaultValidator
from proteovault.database import ProteoVaultDatabase


class VaultPathUpdater:
    """Replace analysis-root paths without changing imported measurements."""

    def __init__(self, database: str | Path, pipeline_output: str | Path,
                 fasta: str | Path | None = None):
        self.database = Path(database).expanduser().resolve()
        self.pipeline_output = Path(pipeline_output).expanduser().resolve()
        self.fasta = Path(fasta).expanduser().resolve() if fasta else None

    def plan(self) -> dict:
        if not self.database.is_file():
            raise FileNotFoundError(self.database)
        if not self.pipeline_output.is_dir():
            raise NotADirectoryError(self.pipeline_output)
        if self.fasta and not self.fasta.is_file():
            raise FileNotFoundError(self.fasta)
        report = VaultValidator(self.database).validate()
        if not report["valid"]:
            raise ValueError("Invalid database: " + "; ".join(report["errors"]))
        with ProteoVaultDatabase(self.database, read_only=True) as vault:
            old_root = vault.connection.execute("SELECT root_path FROM analysis").fetchone()[0]
            sources = vault.connection.execute(
                "SELECT source_id,path,relative_path,format FROM source_files"
            ).fetchall()
        old_root = Path(old_root)
        if not old_root.is_absolute():
            old_root = (self.database.parent / old_root).resolve()
        changes = []
        for source_id, old_path, relative, source_format in sources:
            if source_format == "fasta":
                if self.fasta is not None:
                    new_path = self.fasta
                elif Path(old_path).is_relative_to(old_root):
                    new_path = self.pipeline_output / Path(old_path).relative_to(old_root)
                else:
                    continue
            else:
                relative_path = Path(relative)
                if relative_path.is_absolute() or ".." in relative_path.parts:
                    raise ValueError(f"Invalid source relative path: {relative}")
                new_path = self.pipeline_output / relative_path
            if not new_path.is_file():
                raise FileNotFoundError(new_path)
            if str(new_path) != old_path:
                changes.append((source_id, old_path, str(new_path)))
        return dict(database=str(self.database), old_root=str(old_root),
                    new_root=str(self.pipeline_output), source_changes=changes,
                    fasta=str(self.fasta) if self.fasta else None)

    def update(self) -> dict:
        plan = self.plan()
        descriptor, filename = tempfile.mkstemp(
            prefix=".proteovault-paths-", suffix=".duckdb", dir=self.database.parent
        )
        os.close(descriptor)
        temporary = Path(filename)
        try:
            with ProteoVaultDatabase(self.database, read_only=True):
                shutil.copyfile(self.database, temporary)
            with ProteoVaultDatabase(temporary) as vault:
                connection = vault.connection
                connection.execute("BEGIN TRANSACTION")
                try:
                    connection.execute("UPDATE analysis SET root_path=?", [plan["new_root"]])
                    for source_id, _old_path, new_path in plan["source_changes"]:
                        connection.execute(
                            "UPDATE source_files SET path=? WHERE source_id=?",
                            [new_path, source_id],
                        )
                    old_prefix = plan["old_root"].rstrip(os.sep) + os.sep
                    rows = connection.execute(
                        "SELECT source_id,row_number,source_file FROM qc_metrics "
                        "WHERE starts_with(source_file, ?)", [old_prefix]
                    ).fetchall()
                    for source_id, row_number, old_path in rows:
                        connection.execute(
                            "UPDATE qc_metrics SET source_file=? "
                            "WHERE source_id=? AND row_number=?",
                            [plan["new_root"] + os.sep + old_path[len(old_prefix):],
                             source_id, row_number],
                        )
                    connection.execute("COMMIT")
                except Exception:
                    connection.execute("ROLLBACK")
                    raise
            report = VaultValidator(temporary).validate()
            if not report["valid"]:
                raise ValueError("Path update failed validation: " + "; ".join(report["errors"]))
            os.replace(temporary, self.database)
            return dict(database=str(self.database), updated_sources=len(plan["source_changes"]),
                        old_root=plan["old_root"], new_root=plan["new_root"],
                        updated_qc_paths=len(rows), validation=report)
        finally:
            temporary.unlink(missing_ok=True)
            Path(str(temporary) + ".wal").unlink(missing_ok=True)


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser("update-paths", help="Repair paths after moving pipeline results")
    parser.add_argument("database", type=Path)
    parser.add_argument("--pipeline-output", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, help="New location of the search FASTA")
    parser.add_argument("--dry-run", action="store_true")
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    updater = VaultPathUpdater(args.database, args.pipeline_output, args.fasta)
    print(json.dumps(updater.plan() if args.dry_run else updater.update(), indent=2))
    return 0
