"""Integration tests for the installed command handlers."""

import gzip
import json
from pathlib import Path

import duckdb
import pytest

from proteovault import (
    VaultBuilder,
    VaultFastaAdder,
    VaultInspector,
    VaultLayerAdder,
    VaultValidator,
)
from proteovault.cli.main import main


class TestBuildCLI:
    def test_cli_and_dry_run(self, pipeline, tmp_path, capsys):
        output = tmp_path / "not_created" / "results.duckdb"
        args = ["build", "--pipeline-output", str(pipeline), "--output", str(output)]
        assert main(args + ["--dry-run"]) == 0
        assert not output.parent.exists()
        plan = json.loads(capsys.readouterr().out)
        assert len(plan["artifacts"]) == 4
        assert main(args) == 0
        capsys.readouterr()
        assert main(["validate", str(output)]) == 0
        assert json.loads(capsys.readouterr().out)["valid"]
        assert main(["inspect", str(output)]) == 0
        inspection = capsys.readouterr().out
        assert "Table counts" in inspection
        assert "Sources" in inspection
        assert "psms" in inspection
        assert any(
            line.strip().split() == ["psms", "1"] for line in inspection.splitlines()
        )
        assert main(["inspect", str(output), "--json"]) == 0
        assert json.loads(capsys.readouterr().out)["table_counts"]["psms"] == 1
        assert main(args) == 1
        assert "already exists" in capsys.readouterr().err

    def test_cli_logs_progress_to_stderr_only(self, pipeline, tmp_path, capsys):
        output = tmp_path / "logging.duckdb"
        args = ["build", "--pipeline-output", str(pipeline), "--output", str(output)]
        assert main(args) == 0
        captured = capsys.readouterr()
        assert json.loads(captured.out)["database"] == str(output)
        for message in (
            "Discovering pipeline outputs",
            "Importing file 1/4",
            "Imported ",
            "Validating the completed database",
            "Validation passed",
            "Build completed",
        ):
            assert message in captured.err
        assert captured.err.count("Importing file 1/4") == 1
        assert main(["inspect", str(output)]) == 0
        captured = capsys.readouterr()
        assert "Table counts" in captured.out
        assert "psms" in captured.out
        assert "Inspection complete" in captured.err
        assert main(["validate", str(output)]) == 0
        captured = capsys.readouterr()
        assert json.loads(captured.out)["valid"]
        assert "Validation passed" in captured.err
        assert main(args) == 1
        captured = capsys.readouterr()
        assert not captured.out
        assert "build failed" in captured.err
