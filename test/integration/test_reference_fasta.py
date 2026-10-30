"""Integration tests for linked imports, provenance, and atomic publication."""

import csv
import gzip
import json
import shutil
import xml.etree.ElementTree as ET
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

TEST_FIXTURES = Path(__file__).parents[1] / "fixtures"


class TestReferenceFasta:
    def test_build_imports_search_fasta_for_existing_proteins(self, pipeline, tmp_path):
        fasta = tmp_path / "search.fasta"
        shutil.copyfile(TEST_FIXTURES / "search.fasta", fasta)
        output = tmp_path / "with-fasta.duckdb"
        plan = VaultBuilder(pipeline, output, fasta=fasta).plan()
        assert plan["artifacts"][-1]["format"] == "fasta"
        assert plan["artifacts"][-1]["relative_path"] == "reference_fasta/search.fasta"
        result = VaultBuilder(pipeline, output, fasta=fasta).build()
        assert result["validation"]["valid"]
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute(
                "SELECT accession, sequence FROM proteins ORDER BY accession"
            ).fetchall() == [
                ("P1", "PEPTIDE"),
                ("P2", "MPEPTIDER"),
                ("P3", None),
            ]
            assert connection.execute(
                "SELECT relative_path, path FROM source_files WHERE format='fasta'"
            ).fetchall() == [("reference_fasta/search.fasta", str(fasta))]
            assert connection.execute(
                "SELECT row_count FROM import_counts WHERE entity='protein_entries'"
            ).fetchall() == [(2,)]

    def test_build_accepts_gzipped_fasta_and_rejects_wrong_database(
        self, pipeline, tmp_path
    ):
        fasta = tmp_path / "search.fasta.gz"
        with gzip.open(fasta, "wt") as handle:
            handle.write(">P2\nMPEPTIDER\n")
        output = tmp_path / "gzipped.duckdb"
        assert (
            main(
                [
                    "build",
                    "--pipeline-output",
                    str(pipeline),
                    "--output",
                    str(output),
                    "--fasta",
                    str(fasta),
                ]
            )
            == 0
        )
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute(
                    "SELECT sequence FROM proteins WHERE accession='P2'"
                ).fetchone()[0]
                == "MPEPTIDER"
            )
        wrong = tmp_path / "wrong.fasta"
        wrong.write_text(">UNRELATED\nAAAAAA\n")
        rejected = tmp_path / "rejected.duckdb"
        with pytest.raises(ValueError, match="No FASTA accessions match"):
            VaultBuilder(pipeline, rejected, fasta=wrong).build()
        assert not rejected.exists()

    def test_build_rejects_conflicting_fasta_sequence(self, pipeline, tmp_path):
        fasta = tmp_path / "conflict.fasta"
        fasta.write_text(">P1\nWRONG\n")
        output = tmp_path / "conflict.duckdb"
        with pytest.raises(ValueError, match="conflicts with existing sequence"):
            VaultBuilder(pipeline, output, fasta=fasta).build()
        assert not output.exists()

    def test_add_fasta_fills_missing_sequences_and_is_idempotent(
        self, pipeline, tmp_path
    ):
        output = tmp_path / "existing.duckdb"
        VaultBuilder(pipeline, output).build()
        fasta = tmp_path / "search.fasta"
        fasta.write_text(">P1\nPEPTIDE\n" ">P2\nMPEPTIDER\n" ">UNUSED\nAAAAAA\n")
        adder = VaultFastaAdder(database=output, fasta=fasta)
        assert adder.plan()["proteins_missing_sequence"] == 2
        result = adder.add()
        assert result["imported_files"] == 1
        assert result["matched_proteins"] == 2
        assert result["new_sequences"] == 1
        assert result["proteins_missing_sequence"] == 1
        assert result["validation"]["valid"]
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute(
                "SELECT accession, sequence FROM proteins ORDER BY accession"
            ).fetchall() == [
                ("P1", "PEPTIDE"),
                ("P2", "MPEPTIDER"),
                ("P3", None),
            ]
            assert connection.execute(
                "SELECT relative_path FROM source_files WHERE format='fasta'"
            ).fetchall() == [("reference_fasta/search.fasta",)]
        assert "fasta" in VaultInspector(output).summary()["layers"]
        before = output.read_bytes()
        assert adder.add()["status"] == "already_imported"
        assert output.read_bytes() == before

    def test_add_fasta_rejects_mismatch_without_changing_vault(
        self, pipeline, tmp_path
    ):
        output = tmp_path / "existing.duckdb"
        VaultBuilder(pipeline, output).build()
        wrong = tmp_path / "wrong.fasta"
        wrong.write_text(">P1\nWRONG\n")
        before = output.read_bytes()
        with pytest.raises(ValueError, match="conflicts with existing sequence"):
            VaultFastaAdder(database=output, fasta=wrong).add()
        assert output.read_bytes() == before
