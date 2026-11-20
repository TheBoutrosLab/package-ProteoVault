"""Build a portable analysis database from published DIA pipeline artifacts."""

from __future__ import annotations

from typing import TYPE_CHECKING
import json
import logging
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import duckdb

from proteovault.cli.inspect import VaultInspector
from proteovault.cli.validate import VaultValidator
from proteovault.database import ProteoVaultDatabase
from proteovault.discovery import PipelineOutputDiscovery
from proteovault.importers.coordinator import ImportCoordinator
from proteovault.schema import SCHEMA_VERSION
from proteovault.version import __version__

if TYPE_CHECKING:
    from argparse import Namespace, ArgumentParser, _SubParsersAction

logger = logging.getLogger(__name__)


class VaultBuilder:
    """Discover inputs, coordinate format importers, and publish a validated database."""

    def __init__(
        self,
        pipeline_output: str | Path,
        output: str | Path,
        force: bool = False,
        fasta: str | Path | None = None,
        uniprot_idmapping: str | Path | None = None,
        tmt_channel_csv: str | Path | None = None,
        workflow_type: str | None = None,
    ):
        self.pipeline_output = Path(pipeline_output).expanduser().resolve()
        self.output = Path(output).expanduser().absolute()
        self.force = force
        self.fasta = Path(fasta).expanduser().resolve() if fasta else None
        self.uniprot_idmapping = (
            Path(uniprot_idmapping).expanduser().resolve() if uniprot_idmapping else None
        )
        self.tmt_channel_csv = (
            Path(tmt_channel_csv).expanduser().resolve() if tmt_channel_csv else None
        )
        self.workflow_type = workflow_type

    def plan(self) -> dict:
        if self.output.exists() and not self.force:
            raise FileExistsError(
                f"Output already exists: {self.output}; use --force to replace it"
            )
        logger.info("Discovering pipeline outputs in %s", self.pipeline_output)
        plan = PipelineOutputDiscovery(self.pipeline_output, self.workflow_type).discover()
        if plan["workflow_type"] == "DDA" and plan["acquisition_mode"] == "TMT":
            if self.tmt_channel_csv is None or not self.tmt_channel_csv.is_file():
                raise ValueError("DDA TMT builds require --tmt-channel-csv")
            plan["artifacts"].insert(0, dict(
                path=str(self.tmt_channel_csv), format="tmt_channel_manifest",
                stage="design", sample_name=None,
                relative_path=f"tmt_channel_manifest/{self.tmt_channel_csv.name}",
            ))
        elif self.tmt_channel_csv is not None:
            raise ValueError("--tmt-channel-csv is only valid for DDA TMT builds")
        if self.fasta is not None:
            if not self.fasta.is_file():
                raise FileNotFoundError(f"FASTA file does not exist: {self.fasta}")
            plan["artifacts"].append(
                dict(
                    path=str(self.fasta),
                    format="fasta",
                    stage="reference",
                    sample_name=None,
                    relative_path=f"reference_fasta/{self.fasta.name}",
                )
            )
        if self.uniprot_idmapping is not None:
            if not self.uniprot_idmapping.is_file():
                raise FileNotFoundError(
                    f"UniProt ID mapping file does not exist: {self.uniprot_idmapping}"
                )
            plan["artifacts"].append(
                dict(
                    path=str(self.uniprot_idmapping),
                    format="uniprot_idmapping",
                    stage="annotation",
                    sample_name=None,
                    relative_path=f"uniprot_idmapping/{self.uniprot_idmapping.name}",
                )
            )
        logger.info(
            "Discovered %d files for %s (%s %s)",
            len(plan["artifacts"]),
            plan["dataset_id"],
            plan["pipeline"],
            plan["pipeline_version"],
        )
        for warning in plan["warnings"]:
            logger.warning("%s", warning)
        for artifact in plan["artifacts"]:
            if Path(artifact["path"]).resolve() == self.output.resolve():
                raise ValueError("Output would replace a source artifact")
        return {**plan, "output": str(self.output)}

    def build(self) -> dict:
        started = perf_counter()
        plan = self.plan()
        logger.info("Building database for %s at %s", plan["dataset_id"], self.output)
        self.output.parent.mkdir(parents=True, exist_ok=True)
        descriptor, name = tempfile.mkstemp(
            prefix=".proteovault-", suffix=".duckdb", dir=self.output.parent
        )
        os.close(descriptor)
        temporary = Path(name)
        temporary.unlink()  # DuckDB initializes a new file rather than opening an empty file.
        try:
            with ProteoVaultDatabase(temporary) as database:
                database.initialize()
                database.connection.execute("BEGIN TRANSACTION")
                try:
                    for key, value in {
                        "format": "ProteoVault",
                        "schema_version": SCHEMA_VERSION,
                        "workflow_type": plan["workflow_type"],
                        "producer_version": __version__,
                        "duckdb_version": duckdb.__version__,
                        "created_at": datetime.now(timezone.utc).isoformat(),
                        "warnings": json.dumps(plan["warnings"]),
                    }.items():
                        database.add("metadata", key=key, value=value)
                    database.add(
                        "analysis",
                        analysis_id=1,
                        dataset_id=plan["dataset_id"],
                        pipeline=plan["pipeline"],
                        pipeline_version=plan["pipeline_version"],
                        root_path=plan["root"],
                    )
                    import_warnings = ImportCoordinator(
                        database, self.pipeline_output, plan["workflow_type"],
                        plan.get("acquisition_mode"),
                    ).import_artifacts(plan["artifacts"])
                    if import_warnings:
                        plan["warnings"].extend(import_warnings)
                        database.connection.execute(
                            "UPDATE metadata SET value=? WHERE key='warnings'",
                            [json.dumps(plan["warnings"])],
                        )
                    database.connection.execute("COMMIT")
                except Exception:
                    database.connection.execute("ROLLBACK")
                    raise
            logger.info("Validating the completed database before publication")
            report = VaultValidator(temporary).validate()
            if not report["valid"]:
                raise ValueError(
                    "Built database failed validation: " + "; ".join(report["errors"])
                )
            if self.force:
                os.replace(temporary, self.output)
            else:
                # Atomic creation without overwriting a concurrently created destination.
                os.link(temporary, self.output)
                temporary.unlink()
            counts = VaultInspector(self.output).table_counts()
            logger.info(
                "Build completed in %.2f s: %s (%d runs, %d PSMs, %d protein measurements)",
                perf_counter() - started,
                self.output,
                counts["runs"],
                counts["psms"],
                counts["protein_quant"],
            )
            return {
                "database": str(self.output),
                "validation": report,
                "table_counts": counts,
            }
        finally:
            temporary.unlink(missing_ok=True)
            Path(str(temporary) + ".wal").unlink(missing_ok=True)

    @staticmethod
    def importer_class(format_name: str):
        return ImportCoordinator.importer_class(format_name)


def add_parser(subparsers: _SubParsersAction[ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "build", help="Import core quantification, validated PSMs, and protein groups"
    )
    parser.add_argument("--pipeline-output", required=True, type=Path)
    parser.add_argument("--workflow-type", choices=("DIA", "DDA"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--fasta",
        type=Path,
        help="Target-decoy FASTA used by the search; import sequences for identified proteins",
    )
    parser.add_argument(
        "--uniprot-idmapping",
        type=Path,
        help="Local UniProt idmapping.gz archive for existing UniProt proteins",
    )
    parser.add_argument(
        "--tmt-channel-csv", type=Path,
        help="TMT channel.csv with plex_id,channel,description columns",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace output only after a successful build",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show input artifacts without writing a database",
    )
    parser.set_defaults(handler=run)


def run(args: Namespace) -> int:
    builder = VaultBuilder(
        args.pipeline_output,
        args.output,
        force=args.force,
        fasta=args.fasta,
        uniprot_idmapping=args.uniprot_idmapping,
        tmt_channel_csv=args.tmt_channel_csv,
        workflow_type=args.workflow_type,
    )
    result = builder.plan() if args.dry_run else builder.build()
    if args.dry_run:
        logger.info("Dry run complete; no database was written")
    print(json.dumps(result, indent=2))
    return 0
