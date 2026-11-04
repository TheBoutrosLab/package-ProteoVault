"""Add UniProt cross-references to an existing ProteoVault database."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import os
import shutil
import tempfile
from pathlib import Path

from proteovault.cli.validate import VaultValidator
from proteovault.database import ProteoVaultDatabase
from proteovault.importers.uniprot_idmapping import UniProtIdMappingImporter
from proteovault.manifest import SourceManifest

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction


class VaultUniProtMappingAdder:
    """Atomically add one local mapping archive to a validated vault."""

    def __init__(self, database: str | Path, mapping: str | Path):
        self.database = Path(database).expanduser().resolve()
        self.mapping = Path(mapping).expanduser().resolve()

    def plan(self) -> dict:
        if not self.database.is_file():
            raise FileNotFoundError(f"Database does not exist: {self.database}")
        if not self.mapping.is_file():
            raise FileNotFoundError(f"UniProt mapping file does not exist: {self.mapping}")
        if self.database == self.mapping:
            raise ValueError("Mapping archive must differ from the database")
        report = VaultValidator(self.database).validate()
        if not report["valid"]:
            raise ValueError("Cannot extend an invalid database: " + "; ".join(report["errors"]))
        relative = f"uniprot_idmapping/{self.mapping.name}"
        stat = self.mapping.stat()
        with ProteoVaultDatabase(self.database, read_only=True) as vault:
            existing = vault.connection.execute(
                "SELECT path, relative_path, size_bytes, mtime_ns FROM source_files "
                "WHERE relative_path=? OR path=?",
                [relative, str(self.mapping)],
            ).fetchall()
        if existing:
            if len(existing) != 1 or existing[0] != (
                str(self.mapping), relative, stat.st_size, stat.st_mtime_ns
            ):
                raise ValueError("A different mapping source is already recorded; rebuild the vault")
            status = "already_imported"
        else:
            status = "ready"
        return dict(database=str(self.database), mapping=str(self.mapping),
                    relative_path=relative, size_bytes=stat.st_size,
                    mtime_ns=stat.st_mtime_ns, status=status)

    def add(self) -> dict:
        plan = self.plan()
        if plan["status"] == "already_imported":
            return dict(database=str(self.database), mapping=str(self.mapping),
                        imported_files=0, status="already_imported")
        descriptor, filename = tempfile.mkstemp(
            prefix=".proteovault-idmapping-", suffix=".duckdb", dir=self.database.parent
        )
        os.close(descriptor)
        temporary = Path(filename)
        database_stat = self.database.stat()
        signature = (database_stat.st_ino, database_stat.st_size, database_stat.st_mtime_ns)
        try:
            with ProteoVaultDatabase(self.database, read_only=True):
                shutil.copyfile(self.database, temporary)
            with ProteoVaultDatabase(temporary) as vault:
                vault.load_layout()
                vault.connection.execute("BEGIN TRANSACTION")
                try:
                    artifact = dict(
                        path=str(self.mapping), format="uniprot_idmapping",
                        stage="annotation", sample_name=None,
                        relative_path=plan["relative_path"],
                    )
                    source_id = SourceManifest(self.database.parent).add_source(vault, artifact)
                    counts = UniProtIdMappingImporter(
                        path=self.mapping, database=vault, source_id=source_id,
                        stage="annotation"
                    ).import_data()
                    vault.connection.execute("COMMIT")
                except Exception:
                    vault.connection.execute("ROLLBACK")
                    raise
            report = VaultValidator(temporary).validate()
            if not report["valid"]:
                raise ValueError("Mapping import failed validation: " + "; ".join(report["errors"]))
            current = self.database.stat()
            if (current.st_ino, current.st_size, current.st_mtime_ns) != signature:
                raise ValueError("Database changed during mapping import; no result was published")
            mapping_stat = self.mapping.stat()
            if (mapping_stat.st_size, mapping_stat.st_mtime_ns) != (
                plan["size_bytes"], plan["mtime_ns"]
            ):
                raise ValueError("Mapping archive changed during import; no result was published")
            os.replace(temporary, self.database)
            return dict(database=str(self.database), mapping=str(self.mapping),
                        imported_files=1, matched_proteins=counts["matched_proteins"],
                        unmatched_proteins=counts["unmatched_proteins"],
                        skipped_proteins=counts["skipped_proteins"],
                        identifiers=counts["protein_identifiers"], validation=report)
        finally:
            temporary.unlink(missing_ok=True)
            Path(str(temporary) + ".wal").unlink(missing_ok=True)


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "add-uniprot-idmapping", help="Add local UniProt protein cross-references"
    )
    parser.add_argument("database", type=Path)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    adder = VaultUniProtMappingAdder(args.database, args.mapping)
    print(json.dumps(adder.plan() if args.dry_run else adder.add(), indent=2))
    return 0
