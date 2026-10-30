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


class TestAddFastaCLI:
    def test_add_fasta_cli_and_failed_import_preserve_vault(
        self, pipeline, tmp_path, capsys
    ):
        output = tmp_path / "existing.duckdb"
        VaultBuilder(pipeline, output).build()
        fasta = tmp_path / "search.fasta.gz"
        with gzip.open(fasta, "wt") as handle:
            handle.write(">P2\nMPEPTIDER\n")
        command = ["add-fasta", str(output), "--fasta", str(fasta)]
        assert main(command + ["--dry-run"]) == 0
        assert json.loads(capsys.readouterr().out)["status"] == "ready"
        assert main(command) == 0
        result = json.loads(capsys.readouterr().out)
        assert result["new_sequences"] == 1
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute(
                    "SELECT sequence FROM proteins WHERE accession='P2'"
                ).fetchone()[0]
                == "MPEPTIDER"
            )
        wrong = tmp_path / "wrong.fasta"
        wrong.write_text(">UNRELATED\nAAAAAA\n")
        before = output.read_bytes()
        with pytest.raises(ValueError, match="different reference FASTA"):
            VaultFastaAdder(database=output, fasta=wrong).add()
        assert output.read_bytes() == before
