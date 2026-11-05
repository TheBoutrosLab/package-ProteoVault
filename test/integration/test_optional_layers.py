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


class TestOptionalLayers:
    def test_pin_features_are_preserved_without_guessed_links(self, pipeline, tmp_path):
        workflow = next(pipeline.glob("*/*/main_workflow"))
        path = (
            workflow
            / "validate-ProteomicPSM-2.2.0/S1/MSBooster/output/engine_MSFragger/run_1-fraction_1_edited.pin"
        )
        path.parent.mkdir(parents=True)
        shutil.copyfile(TEST_FIXTURES / "msbooster.pin", path)
        output = tmp_path / "pin.duckdb"
        VaultBuilder(pipeline, output).build()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute("SELECT count(*) FROM pin_features").fetchone()[0]
                == 0
            )
        VaultLayerAdder(output, ["msbooster"]).add()
        with duckdb.connect(str(output), read_only=True) as connection:
            assert connection.execute(
                "SELECT spectrum_reference, feature_name, value FROM pin_features"
            ).fetchall() == [("S1.7.7.2", "deltaRTLOESS", 0.2)]

    def test_probability_validation_and_zero(self, pipeline, tmp_path):
        report = next(pipeline.rglob("*-report.tsv"))
        report.write_text(report.read_text().replace("0.005", "1.1"))
        output = tmp_path / "invalid.duckdb"
        with pytest.raises(ValueError, match="invalid confidence"):
            VaultBuilder(pipeline, output).build()
        assert not output.exists()

    def test_optional_sources_are_not_read_by_core_build(self, pipeline, tmp_path):
        source = next(pipeline.rglob("S1_MSFragger.idXML"))
        source.write_text("not XML")
        output = tmp_path / "core.duckdb"
        assert VaultBuilder(pipeline, output).build()["validation"]["valid"]
        assert VaultInspector(output).summary()["layers"] == ["core"]
        before = output.read_bytes()
        with pytest.raises(ET.ParseError):
            VaultLayerAdder(output, ["raw-search"]).add()
        assert output.read_bytes() == before
        assert not list(tmp_path.glob(".proteovault-layer-*"))

    def test_layer_rejects_changed_imported_source(self, pipeline, tmp_path):
        output = tmp_path / "changed.duckdb"
        VaultBuilder(pipeline, output).build()
        VaultLayerAdder(output, ["raw-search"]).add()
        source = next(pipeline.rglob("S1_MSFragger.idXML"))
        source.write_text(source.read_text() + "\n")
        before = output.read_bytes()
        with pytest.raises(ValueError, match="Previously imported source changed"):
            VaultLayerAdder(output, ["raw-search"]).add()
        assert output.read_bytes() == before

    def test_layer_rejects_other_dataset(self, pipeline, tmp_path):
        import shutil

        other = tmp_path / "other"
        shutil.copytree(pipeline, other)
        version = next(other.glob("metapipeline-ProteomicDIA-*"))
        (version / "TEST").rename(version / "OTHER")
        output = tmp_path / "same.duckdb"
        VaultBuilder(pipeline, output).build()
        with pytest.raises(ValueError, match="different dataset"):
            VaultLayerAdder(output, ["raw-search"], pipeline_output=other).add()

    def test_layer_root_alias_preserves_source_identity(self, pipeline, tmp_path):
        output = tmp_path / "alias.duckdb"
        VaultBuilder(pipeline, output).build()
        version = next(pipeline.glob("metapipeline-ProteomicDIA-*"))
        VaultLayerAdder(output, ["raw-search"], pipeline_output=version).add()
        assert VaultLayerAdder(output, ["raw-search"]).add()["imported_files"] == 0

    def test_multiple_layers_publish_together(self, pipeline, tmp_path):
        workflow = next(pipeline.glob("*/*/main_workflow"))
        pin = (
            workflow
            / "validate-ProteomicPSM-2.2.0/S1/MSBooster/output/engine_MSFragger/run_1-fraction_1_edited.pin"
        )
        pin.parent.mkdir(parents=True)
        pin.write_text("broken PIN header\n")
        output = tmp_path / "combined.duckdb"
        VaultBuilder(pipeline, output).build()
        before = output.read_bytes()
        with pytest.raises(ValueError, match="PIN requires"):
            VaultLayerAdder(output, ["raw-search", "msbooster"]).add()
        assert output.read_bytes() == before
        pin.write_text(
            "SpecId\tLabel\tScanNr\tdeltaRTLOESS\tPeptide\tProteins\nS1.7.7.2\t1\t7\t0.2\tK.PEPTIDE.R\tP1\n"
        )
        result = VaultLayerAdder(output, ["raw-search", "msbooster"]).add()
        assert result["imported_files"] == 2
        assert result["table_counts"]["psms"] == 2
        assert result["table_counts"]["pin_features"] == 1
        assert VaultInspector(output).summary()["layers"] == [
            "core",
            "msbooster",
            "raw-search",
        ]

    def test_compact_schema_and_report_quantities_in_core(self, pipeline, tmp_path):
        output = tmp_path / "compact.duckdb"
        built = VaultBuilder(pipeline, output).build()
        assert built["validation"]["valid"]
        assert built["table_counts"]["precursor_quant"] == 4
        with duckdb.connect(str(output), read_only=True) as connection:
            assert (
                connection.execute(
                    "SELECT value FROM metadata WHERE key='schema_version'"
                ).fetchone()[0]
                == "1.1.0"
            )
            assert (
                connection.execute(
                    "SELECT data_type FROM information_schema.columns "
                    "WHERE table_name='psms' AND column_name='psm_id'"
                ).fetchone()[0]
                == "BIGINT"
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM duckdb_constraints() WHERE constraint_type='FOREIGN KEY'"
                ).fetchone()[0]
                == 0
            )
            assert (
                connection.execute(
                    "SELECT count(*) FROM precursor_quant_readable WHERE method='matrix'"
                ).fetchone()[0]
                == 2
            )
        assert VaultInspector(output).summary()["table_counts"]["precursor_quant"] == 4
        assert VaultInspector(output).summary()["layers"] == ["core"]
