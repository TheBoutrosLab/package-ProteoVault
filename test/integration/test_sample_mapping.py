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
