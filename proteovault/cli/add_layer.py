"""Extend an existing vault with explicitly selected optional source layers."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import logging
import os
import re
import shutil
import tempfile
from pathlib import Path
from time import perf_counter

from proteovault.cli.inspect import VaultInspector
from proteovault.cli.validate import VaultValidator
from proteovault.database import ProteoVaultDatabase
from proteovault.discovery import PipelineOutputDiscovery
from proteovault.importers.coordinator import ImportCoordinator

if TYPE_CHECKING:
    from argparse import ArgumentParser, Namespace, _SubParsersAction

logger = logging.getLogger(__name__)


class VaultLayerAdder:
    """Add sources on a temporary copy and publish only a validated result."""

    def __init__(
        self,
        database: str | Path,
        layers: list[str],
        pipeline_output: str | Path | None = None,
    ) -> None:
        self.database = Path(database).expanduser().absolute()
        self.layers = list(dict.fromkeys(layers))
        self.pipeline_output = (
            Path(pipeline_output).expanduser().resolve() if pipeline_output else None
        )
        if not self.layers or any(
            layer not in ("quant-details", "raw-search", "msbooster", "qc")
            for layer in self.layers
        ):
            raise ValueError(
                "Choose at least one optional layer: quant-details, raw-search, msbooster, qc"
            )

    def plan(self) -> dict:
        logger.info("Planning layers %s for %s", ", ".join(self.layers), self.database)
        if not self.database.is_file():
            raise FileNotFoundError(f"Database does not exist: {self.database}")
        report = VaultValidator(self.database).validate()
        if not report["valid"]:
            raise ValueError(
                "Cannot extend an invalid database: " + "; ".join(report["errors"])
            )
        with ProteoVaultDatabase(self.database, read_only=True) as database:
            row = database.connection.execute(
                "SELECT dataset_id,pipeline,pipeline_version,root_path FROM analysis"
            ).fetchone()
            known = {
                relative: (size, mtime, format_name, stage)
                for relative, size, mtime, format_name, stage in database.connection.execute(
                    "SELECT relative_path,size_bytes,mtime_ns,format,stage FROM source_files"
                ).fetchall()
            }
        root = self.pipeline_output or Path(row[3])
        discovery = PipelineOutputDiscovery(root)
        workflow, _ = discovery.analysis()
        source_root = (
            workflow.parents[1]
            if re.fullmatch(r"metapipeline-ProteomicDIA-[0-9].*", Path(row[3]).name)
            else workflow.parents[2]
        )
        artifacts = []
        skipped = []
        for layer in self.layers:
            discovered = discovery.discover(layer=layer)
            if (
                tuple(
                    discovered[key]
                    for key in ("dataset_id", "pipeline", "pipeline_version")
                )
                != row[:3]
            ):
                raise ValueError(
                    "Optional layer belongs to a different dataset or pipeline version"
                )
            for artifact in discovered["artifacts"]:
                path = Path(artifact["path"])
                relative = path.relative_to(source_root).as_posix()
                stat = path.stat()
                if relative in known:
                    expected = (
                        stat.st_size,
                        stat.st_mtime_ns,
                        artifact["format"],
                        artifact["stage"],
                    )
                    if known[relative] != expected:
                        raise ValueError(
                            f"Previously imported source changed: {path}; rebuild the database"
                        )
                    skipped.append(artifact)
                else:
                    artifacts.append(artifact)
        logger.info(
            "Layer plan: %d files to import, %d already present",
            len(artifacts),
            len(skipped),
        )
        return dict(
            database=str(self.database),
            root=str(source_root),
            layers=self.layers,
            artifacts=artifacts,
            skipped=skipped,
        )

    def add(self) -> dict:
        started = perf_counter()
        plan = self.plan()
        if not plan["artifacts"]:
            logger.info("All requested layer files are already imported")
            return dict(
                database=str(self.database),
                imported_files=0,
                skipped_files=len(plan["skipped"]),
                layers=self.layers,
            )
        descriptor, name = tempfile.mkstemp(
            prefix=".proteovault-layer-", suffix=".duckdb", dir=self.database.parent
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
            logger.info("Copying database before adding optional layers")
            with ProteoVaultDatabase(self.database, read_only=True):
                shutil.copyfile(self.database, temporary)
            with ProteoVaultDatabase(temporary) as database:
                database.load_layout()
                database.connection.execute("BEGIN TRANSACTION")
                try:
                    ImportCoordinator(database, plan["root"]).import_artifacts(
                        plan["artifacts"]
                    )
                    database.connection.execute("COMMIT")
                except Exception:
                    database.connection.execute("ROLLBACK")
                    raise
            report = VaultValidator(temporary).validate()
            if not report["valid"]:
                raise ValueError(
                    "Layer import failed validation: " + "; ".join(report["errors"])
                )
            current = self.database.stat()
            if (current.st_ino, current.st_size, current.st_mtime_ns) != signature:
                raise ValueError(
                    "Database changed during layer import; no result was published"
                )
            os.replace(temporary, self.database)
            logger.info(
                "Added %d files in %.2f s",
                len(plan["artifacts"]),
                perf_counter() - started,
            )
            return dict(
                database=str(self.database),
                imported_files=len(plan["artifacts"]),
                skipped_files=len(plan["skipped"]),
                layers=self.layers,
                validation=report,
                table_counts=VaultInspector(self.database).table_counts(),
            )
        finally:
            temporary.unlink(missing_ok=True)
            Path(str(temporary) + ".wal").unlink(missing_ok=True)


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "add-layer",
        help="Add optional quantities, search results, QC metrics, or MSBooster features",
    )
    parser.add_argument("database", type=Path)
    parser.add_argument(
        "--layer",
        choices=("quant-details", "raw-search", "msbooster", "qc"),
        action="append",
        required=True,
        help="Optional layer to import; repeat to select multiple layers",
    )
    parser.add_argument(
        "--pipeline-output",
        type=Path,
        help="Override the recorded pipeline output root after relocation",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    adder = VaultLayerAdder(
        args.database, args.layer, pipeline_output=args.pipeline_output
    )
    print(json.dumps(adder.plan() if args.dry_run else adder.add(), indent=2))
    return 0
