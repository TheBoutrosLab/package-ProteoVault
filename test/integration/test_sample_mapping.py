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


class TestSampleMapping:
    def test_published_quantitation_manifest_sets_canonical_samples(
        self, pipeline, tmp_path
    ):
        workflow = pipeline / "metapipeline-ProteomicDIA-1.3.0/TEST/main_workflow"
        manifests = workflow / "manifests"
        manifests.mkdir()
        shutil.copyfile(TEST_FIXTURES / "quantitation.csv", manifests / "quantitation.csv")
        report = next(pipeline.rglob("*-report.tsv"))
        with report.open("a") as handle:
            handle.write("P1;P2\tPEPTIDE\tPEPTIDE\t2\tS2\t29\t30\t0.005\t0.002\n")
        qc = workflow / "qc-ProteomicDIA-1.0.0/TEST/output/TEST_QC_metrics.tsv"
        qc.parent.mkdir(parents=True)
        qc.write_text(
            (TEST_FIXTURES / "cohort_qc.tsv").read_text().replace("\tS1\t", "\tS2\t")
        )
        output = tmp_path / "mapped.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute(
                "SELECT sample_name FROM samples ORDER BY sample_name"
            ).fetchall() == [("Biological_sample_2",), ("S1",)]
            assert connection.execute(
                "SELECT r.run_name, s.sample_name FROM runs r "
                "JOIN samples s USING(sample_id) ORDER BY r.run_name"
            ).fetchall() == [("S1", "S1"), ("S2", "Biological_sample_2")]
            assert connection.execute(
                "SELECT sample_name FROM protein_quant_long " "WHERE quantity=200"
            ).fetchall() == [("Biological_sample_2",)]
            assert connection.execute(
                "SELECT DISTINCT sample_name FROM qc_metrics WHERE sample_name != ''"
            ).fetchall() == [("Biological_sample_2",)]
            assert connection.execute(
                "SELECT count(*) FROM qc_metrics WHERE sample_name = ''"
            ).fetchone()[0] == 1
        VaultLayerAdder(output, ["quant-details"]).add()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute(
                "SELECT sample_name FROM precursor_quant_long "
                "WHERE method='Precursor.Quantity' AND quantity=29"
            ).fetchall() == [("Biological_sample_2",)]

    def test_legacy_mzml_manifest_remains_supported(self, tmp_path):
        manifest = tmp_path / "quantitation.csv"
        manifest.write_text("sample_name,mzML\nShort,/raw/long.mzML\n")
        from proteovault.sample_manifest import SampleManifest

        assert SampleManifest(manifest).sample_name("long") == "Short"

    def test_qc_name_without_manifest_must_match_quantification(self, pipeline, tmp_path):
        workflow = pipeline / "metapipeline-ProteomicDIA-1.3.0/TEST/main_workflow"
        qc = workflow / "qc-ProteomicDIA-1.0.0/TEST/output/TEST_QC_metrics.tsv"
        qc.parent.mkdir(parents=True)
        qc.write_text(
            (TEST_FIXTURES / "cohort_qc.tsv").read_text().replace(
                "\tS1\t", "\tShortSample\t"
            )
        )
        output = tmp_path / "unmapped.duckdb"
        with pytest.raises(ValueError, match="does not match a quantified sample"):
            VaultBuilder(pipeline, output).build()
        assert not output.exists()

    def test_ambiguous_manifest_stem_is_rejected(self, pipeline, tmp_path):
        workflow = pipeline / "metapipeline-ProteomicDIA-1.3.0/TEST/main_workflow"
        manifests = workflow / "manifests"
        manifests.mkdir()
        (manifests / "quantitation.csv").write_text(
            "sample_name,mzML\nA,/raw/S1.mzML\nB,/other/S1.mzML\n"
        )
        from proteovault.sample_manifest import SampleManifest

        with pytest.raises(ValueError, match="Ambiguous sample mapping"):
            SampleManifest(manifests / "quantitation.csv").sample_name("S1")
