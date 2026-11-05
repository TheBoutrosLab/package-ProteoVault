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


class TestCoreBuild:
    def test_identified_peptide_can_lack_quantitation(self, pipeline, tmp_path):
        matrix = next(pipeline.rglob("*-report.pr_matrix.tsv"))
        matrix.write_text(
            matrix.read_text().replace(
                "P1;P2\tPEPTIDE\tPEPTIDE\t2\t20\t",
                "P1;P2\tPEPTIDE\tPEPTIDE\t2\t\t",
            )
        )
        report = next(pipeline.rglob("*-report.tsv"))
        report.write_text(report.read_text().replace(
            "\tS1\t19\t20\t", "\tS1\t\t\t"
        ))
        output = tmp_path / "unquantified.duckdb"
        assert VaultBuilder(pipeline, output).build()["validation"]["valid"]
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute("SELECT count(*) FROM peptide_proteins").fetchone()[0]
                == 2
            )
            assert (
                connection.execute("SELECT count(*) FROM psm_proteins").fetchone()[0]
                == 2
            )
            assert (
                connection.execute("SELECT count(*) FROM precursor_groups").fetchone()[0]
                == 0
            )

    def test_quantitation_does_not_create_identification_evidence(
        self, pipeline, tmp_path
    ):
        identification = next(pipeline.rglob("*_filter.idXML"))
        identification.unlink()
        output = tmp_path / "quantitation_only.duckdb"
        assert VaultBuilder(pipeline, output).build()["validation"]["valid"]
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute("SELECT count(*) FROM peptide_proteins").fetchone()[0]
                == 0
            )
            assert (
                connection.execute("SELECT count(*) FROM precursor_groups").fetchone()[0]
                == 2
            )

    def test_build_relations_and_nulls(self, pipeline, tmp_path):
        output = tmp_path / "results.duckdb"
        result = VaultBuilder(pipeline, output).build()
        assert result["validation"]["valid"]
        assert result["table_counts"]["psms"] == 1
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute("SELECT count(*) FROM protein_quant").fetchone()[0]
                == 4
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM protein_quant WHERE quantity IS NULL"
                ).fetchone()[0]
                == 1
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM protein_quant WHERE quantity=0"
                ).fetchone()[0]
                == 1
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM protein_group_members"
                ).fetchone()[0]
                == 11
            )
            assert (
                connection.execute("SELECT count(*) FROM psm_proteins").fetchone()[0]
                == 2
            )
            assert (
                connection.execute(
                    "SELECT value FROM psm_scores_readable WHERE name='q-value'"
                ).fetchone()[0]
                == 0.005
            )
            assert (
                connection.execute("SELECT count(*) FROM psm_features").fetchone()[0]
                == 0
            )
            assert (
                connection.execute("SELECT count(*) FROM psm_metadata").fetchone()[0]
                == 0
            )
            assert (
                connection.execute(
                    "SELECT value FROM psm_scores_readable WHERE name='Posterior Error Probability'"
                ).fetchone()[0]
                == 0.002
            )
            assert (
                connection.execute("SELECT stage FROM psms_readable").fetchone()[0]
                == "validation"
            )
            assert connection.execute("SELECT count(*) FROM runs").fetchone()[0] == 2
        assert len(VaultInspector(output).sources()) == 5

    def test_force_failure_preserves_existing(self, pipeline, tmp_path):
        output = tmp_path / "results.duckdb"
        VaultBuilder(pipeline, output).build()
        before = output.read_bytes()
        source = next(pipeline.rglob("*-report.pg_matrix.tsv"))
        source.write_text("Protein.Group\tS1.mzML\nP1\tbroken\n")
        with pytest.raises(ValueError):
            VaultBuilder(pipeline, output, force=True).build()
        assert output.read_bytes() == before
        assert not list(tmp_path.glob(".proteovault-*"))

    def test_ambiguous_layout_and_bad_database(self, pipeline, tmp_path):
        (pipeline / "metapipeline-ProteomicDIA-9.0.0").mkdir()
        with pytest.raises(ValueError, match="Expected one"):
            VaultBuilder(pipeline, tmp_path / "out.duckdb").plan()
        missing = tmp_path / "missing.duckdb"
        assert not VaultValidator(missing).validate()["valid"]
        assert not missing.exists()
        empty = tmp_path / "empty.duckdb"
        duckdb.connect(str(empty)).close()
        assert not VaultValidator(empty).validate()["valid"]

    def test_repeated_matrix_entity_is_rejected(self, pipeline, tmp_path):
        matrix = next(pipeline.rglob("*-report.pg_matrix.tsv"))
        with matrix.open("a") as handle:
            handle.write("P1;P2\t5\t6\n")
        output = tmp_path / "bad.duckdb"
        with pytest.raises(ValueError, match="repeated matrix entity"):
            VaultBuilder(pipeline, output).build()
        assert not output.exists()

    def test_force_success_and_stable_group_ids(self, pipeline, tmp_path):
        output = tmp_path / "force.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            keys = connection.execute(
                "SELECT group_id FROM protein_groups ORDER BY group_id"
            ).fetchall()
        VaultBuilder(pipeline, output, force=True).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute(
                    "SELECT group_id FROM protein_groups ORDER BY group_id"
                ).fetchall()
                == keys
            )

    def test_diann_cached_entities_survive_batch_flush(self, pipeline, tmp_path):
        report = next(pipeline.rglob("*-report.tsv"))
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        with report.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                [
                    "Protein.Group",
                    "Stripped.Sequence",
                    "Modified.Sequence",
                    "Precursor.Charge",
                    "Run",
                    "Precursor.Quantity",
                    "Precursor.Normalised",
                    "PG.MaxLFQ",
                ]
            )
            for sample in ("S1", "S2"):
                for index in range(2600):
                    sequence = "PEP" + "".join(
                        alphabet[(index // (20**digit)) % 20] for digit in range(3)
                    )
                    writer.writerow(["P1;P2", sequence, sequence, 2, sample, 10, 12, 100])
        output = tmp_path / "batched.duckdb"
        assert VaultBuilder(pipeline, output).build()["validation"]["valid"]
        with duckdb.connect(str(output), read_only=True) as connection:
            count, precursors, total = connection.execute(
                "SELECT count(*),count(DISTINCT precursor_id),sum(quantity) "
                "FROM precursor_quant_readable WHERE method='Precursor.Quantity'"
            ).fetchone()
            assert (count, precursors, total) == (5200, 2600, 52000)
            assert connection.execute(
                "SELECT count(*) FROM protein_quant_readable WHERE method='PG.MaxLFQ'"
            ).fetchone()[0] == 2
