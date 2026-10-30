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


class TestAddLayerCLI:
    def test_add_layer_cli_dry_run_and_idempotency(self, pipeline, tmp_path, capsys):
        output = tmp_path / "layers.duckdb"
        VaultBuilder(pipeline, output).build()
        before = output.read_bytes()
        args = ["add-layer", str(output), "--layer", "raw-search"]
        assert main(args + ["--dry-run"]) == 0
        assert len(json.loads(capsys.readouterr().out)["artifacts"]) == 1
        assert output.read_bytes() == before
        assert main(args) == 0
        assert json.loads(capsys.readouterr().out)["imported_files"] == 1
        after = output.read_bytes()
        assert main(args) == 0
        assert json.loads(capsys.readouterr().out)["imported_files"] == 0
        assert output.read_bytes() == after
        assert VaultInspector(output).summary()["layers"] == ["core", "raw-search"]
        assert VaultValidator(output).validate()["valid"]
