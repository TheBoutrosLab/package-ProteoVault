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


def write_cohort_qc(pipeline):
    workflow = next(pipeline.glob("*/*/main_workflow"))
    output = workflow / "qc-ProteomicDIA-1.0.0/TEST/output/TEST_QC_metrics.tsv"
    output.parent.mkdir(parents=True)
    shutil.copyfile(TEST_FIXTURES / "cohort_qc.tsv", output)
    return output


class TestCohortQC:
    def test_cohort_qc_imported_by_build_when_published(self, pipeline, tmp_path):
        write_cohort_qc(pipeline)
        output = tmp_path / "qc.duckdb"
        result = VaultBuilder(pipeline, output).build()
        assert result["validation"]["valid"]
        assert result["table_counts"]["qc_metrics"] == 3
        assert "qc" in VaultInspector(output).summary()["layers"]
        with duckdb.connect(str(output), read_only=True) as connection:
            rows = connection.execute(
                'SELECT sample_name, metric_name, "group", value, status FROM qc_metrics ORDER BY row_number'
            ).fetchall()
            assert rows == [
                ("S1", "ms2_count", "all", "1200", "available"),
                ("S1", "score_type", "all", "expect", "available"),
                ("", "group_count", "all", None, "unavailable"),
            ]

    def test_cohort_qc_can_be_added_later_and_rejects_wrong_dataset(
        self, pipeline, tmp_path
    ):
        output = tmp_path / "qc-later.duckdb"
        VaultBuilder(pipeline, output).build()
        qc_path = write_cohort_qc(pipeline)
        added = VaultLayerAdder(output, ["qc"]).add()
        assert added["table_counts"]["qc_metrics"] == 3
        assert VaultLayerAdder(output, ["qc"]).add()["imported_files"] == 0
        before = output.read_bytes()
        qc_path.write_text(
            qc_path.read_text().replace("TEST\tacquisition", "OTHER\tacquisition")
        )
        with pytest.raises(ValueError, match="belongs to"):
            VaultBuilder(pipeline, tmp_path / "wrong.duckdb").build()
        assert output.read_bytes() == before
