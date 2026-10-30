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


class TestPSMEvidence:
    def test_stages_preserve_scores(self, pipeline, tmp_path):
        output = tmp_path / "stages.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute("SELECT count(*) FROM psms").fetchone()[0] == 1
        VaultLayerAdder(output, ["raw-search"]).add()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute("SELECT count(*) FROM psms").fetchone()[0] == 2
            assert connection.execute("SELECT count(*) FROM spectra").fetchone()[0] == 1
            assert set(
                connection.execute(
                    "SELECT name,value FROM psm_scores_readable"
                ).fetchall()
            ) == {
                ("hyperscore", 42),
                ("q-value", 0.005),
                ("Posterior Error Probability", 0.002),
            }
            assert (
                connection.execute(
                    "SELECT value FROM psm_features WHERE name='expect'"
                ).fetchone()[0]
                == 0.001
            )

    def test_unsupported_idxml_reference_fails(self, pipeline, tmp_path):
        source = next(pipeline.rglob("S1_MSFragger_filter.idXML"))
        source.write_text(
            source.read_text().replace(
                'protein_refs="PH_0 PH_1"', 'protein_refs="PH_99"'
            )
        )
        with pytest.raises(ValueError, match="Unknown protein reference"):
            VaultBuilder(pipeline, tmp_path / "bad.duckdb").build()

    def test_duplicate_confidence_aliases_collapse(self, pipeline, tmp_path):
        source = next(pipeline.rglob("S1_MSFragger_filter.idXML"))
        document = ET.parse(source)
        hit = document.find(".//PeptideHit")
        ET.SubElement(hit, "UserParam", name="MS:1001493", value="0.002")
        document.write(source, encoding="utf-8", xml_declaration=True)
        output = tmp_path / "aliases.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute("SELECT count(*) FROM psm_scores").fetchone()[0] == 2
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM score_aliases WHERE source_name IN "
                    "('MS:1001493','Posterior Error Probability_score')"
                ).fetchone()[0]
                == 2
            )

    def test_distinct_confidence_values_are_not_merged(self, pipeline, tmp_path):
        source = next(pipeline.rglob("S1_MSFragger_filter.idXML"))
        document = ET.parse(source)
        hit = document.find(".//PeptideHit")
        ET.SubElement(hit, "UserParam", name="MS:1001493", value="0.003")
        document.write(source, encoding="utf-8", xml_declaration=True)
        output = tmp_path / "different-scores.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert set(
                connection.execute(
                    "SELECT name, value FROM psm_scores_readable"
                ).fetchall()
            ) == {
                ("q-value", 0.005),
                ("Posterior Error Probability", 0.002),
                ("MS:1001493", 0.003),
            }
