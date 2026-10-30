"""Add reference sequences to an existing ProteoVault database."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path
from time import perf_counter

from proteovault.cli.validate import VaultValidator
from proteovault.database import ProteoVaultDatabase
from proteovault.importers.fasta import FASTAImporter
from proteovault.manifest import SourceManifest

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

logger = logging.getLogger(__name__)


class VaultFastaAdder:
    """Fill known protein sequences on a validated temporary database copy."""

    def __init__(self, database: str | Path, fasta: str | Path) -> None:
        self.database = Path(database).expanduser().resolve()
        self.fasta = Path(fasta).expanduser().resolve()

    def plan(self) -> dict:
        if not self.database.is_file():
            raise FileNotFoundError(f"Database does not exist: {self.database}")
        if not self.fasta.is_file():
            raise FileNotFoundError(f"FASTA file does not exist: {self.fasta}")
        if self.database == self.fasta:
            raise ValueError("The FASTA path must differ from the database path")
        report = VaultValidator(self.database).validate()
        if not report["valid"]:
            raise ValueError(
                "Cannot extend an invalid database: " + "; ".join(report["errors"])
            )
        relative_path = f"reference_fasta/{self.fasta.name}"
        stat = self.fasta.stat()
        with ProteoVaultDatabase(self.database, read_only=True) as database:
            existing = database.connection.execute(
                "SELECT path, relative_path, size_bytes, mtime_ns "
                "FROM source_files WHERE format='fasta' OR relative_path=?",
                [relative_path],
            ).fetchall()
            missing = database.connection.execute(
                "SELECT count(*) FROM proteins WHERE sequence IS NULL OR sequence=''"
            ).fetchone()[0]
        if existing:
            if len(existing) == 1 and existing[0] == (
                str(self.fasta),
                relative_path,
                stat.st_size,
                stat.st_mtime_ns,
            ):
                status = "already_imported"
            else:
                raise ValueError(
                    "The vault already records a different reference FASTA; "
                    "rebuild it with the intended FASTA"
                )
        else:
            status = "ready"
        return dict(
            database=str(self.database),
            fasta=str(self.fasta),
            relative_path=relative_path,
            fasta_size_bytes=stat.st_size,
            fasta_mtime_ns=stat.st_mtime_ns,
            proteins_missing_sequence=missing,
            status=status,
        )

    def add(self) -> dict:
        started = perf_counter()
        plan = self.plan()
        if plan["status"] == "already_imported":
            logger.info("Reference FASTA is already imported")
            return dict(
                database=str(self.database),
                fasta=str(self.fasta),
                imported_files=0,
                status="already_imported",
            )
        descriptor, name = tempfile.mkstemp(
            prefix=".proteovault-fasta-",
            suffix=".duckdb",
            dir=self.database.parent,
        )
        os.close(descriptor)
        temporary = Path(name)
        original_stat = self.database.stat()
        signature = (
            original_stat.st_ino,
            original_stat.st_size,
            original_stat.st_mtime_ns,
        )
        try:
            logger.info("Copying vault before importing %s", self.fasta)
            with ProteoVaultDatabase(self.database, read_only=True):
                shutil.copyfile(self.database, temporary)
            with ProteoVaultDatabase(temporary) as database:
                database.load_layout()
                database.connection.execute("BEGIN TRANSACTION")
                try:
                    artifact = dict(
                        path=str(self.fasta),
                        format="fasta",
                        stage="reference",
                        sample_name=None,
                        relative_path=plan["relative_path"],
                    )
                    source_id = SourceManifest(self.database.parent).add_source(
                        database,
                        artifact,
                    )
                    importer = FASTAImporter(
                        path=self.fasta,
                        database=database,
                        source_id=source_id,
                        stage="reference",
                        only_existing=True,
                    )
                    counts = importer.import_data()
                    proteins_missing_after = database.connection.execute(
                        "SELECT count(*) FROM proteins "
                        "WHERE sequence IS NULL OR sequence=''"
                    ).fetchone()[0]
                    database.connection.execute("COMMIT")
                except Exception:
                    database.connection.execute("ROLLBACK")
                    raise
            report = VaultValidator(temporary).validate()
            if not report["valid"]:
                raise ValueError(
                    "FASTA import failed validation: " + "; ".join(report["errors"])
                )
            current = self.database.stat()
            if (current.st_ino, current.st_size, current.st_mtime_ns) != signature:
                raise ValueError(
                    "Database changed during FASTA import; no result was published"
                )
            fasta_stat = self.fasta.stat()
            if (fasta_stat.st_size, fasta_stat.st_mtime_ns) != (
                plan["fasta_size_bytes"],
                plan["fasta_mtime_ns"],
            ):
                raise ValueError("FASTA changed during import; no result was published")
            os.replace(temporary, self.database)
            logger.info(
                "Added FASTA sequences for %d proteins in %.2f s",
                counts["protein_entries"],
                perf_counter() - started,
            )
            return dict(
                database=str(self.database),
                fasta=str(self.fasta),
                imported_files=1,
                matched_proteins=counts["protein_entries"],
                new_sequences=plan["proteins_missing_sequence"]
                - proteins_missing_after,
                proteins_missing_sequence=proteins_missing_after,
                validation=report,
            )
        finally:
            temporary.unlink(missing_ok=True)
            Path(str(temporary) + ".wal").unlink(missing_ok=True)


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "add-fasta",
        help="Add search FASTA sequences to proteins in an existing vault",
    )
    parser.add_argument("database", type=Path)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    adder = VaultFastaAdder(database=args.database, fasta=args.fasta)
    print(json.dumps(adder.plan() if args.dry_run else adder.add(), indent=2))
    return 0
