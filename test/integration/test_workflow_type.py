"""Select DIA or DDA from the workflow rather than its directory name."""

import json
import shutil
from pathlib import Path

import duckdb
import pytest

from proteovault import VaultBuilder, VaultLayerAdder
from proteovault.cli.main import main
from test_dda_build import dda_output


def snapshot(path):
    with duckdb.connect(str(path), read_only=True) as connection:
        samples = connection.execute(
            "SELECT sample_name FROM samples ORDER BY sample_name"
        ).fetchall()
        proteins = connection.execute(
            "SELECT sample_name,source_group_id,quantity FROM protein_quant_long "
            "ORDER BY sample_name,source_group_id,quantity"
        ).fetchall()
        qc = connection.execute(
            'SELECT dataset_id,sample_name,metric_name,"group",value '
            'FROM qc_metrics ORDER BY dataset_id,sample_name,metric_name,"group"'
        ).fetchall()
        return samples, proteins, qc


@pytest.mark.parametrize("workflow_type", ["DIA", "DDA"])
def test_renamed_workflow_build_and_later_qc(pipeline, tmp_path, workflow_type, capsys):
    if workflow_type == "DIA":
        root = pipeline
        workflow = next(root.glob("*/*/main_workflow"))
        qc = (workflow / "qc-ProteomicDIA-1.0.0/TEST/output/"
              "TEST_QC_metrics.tsv")
        qc.parent.mkdir(parents=True)
        fixture = Path(__file__).parents[1] / "fixtures/cohort_qc.tsv"
        shutil.copyfile(fixture, qc)
    else:
        root = dda_output(tmp_path, "LFQ")
    original = tmp_path / f"original-{workflow_type}.duckdb"
    assert VaultBuilder(root, original).build()["validation"]["valid"]

    renamed = tmp_path / f"renamed-{workflow_type}"
    shutil.copytree(root, renamed)
    version = next(renamed.iterdir())
    version.rename(renamed / "results")
    qc = next(renamed.rglob("*_QC_metrics.tsv"))

    assert main(["build", "--pipeline-output", str(renamed),
                 "--output", str(tmp_path / "dry-run.duckdb"),
                 "--workflow-type", workflow_type, "--dry-run"]) == 0
    assert json.loads(capsys.readouterr().out)["workflow_type"] == workflow_type

    explicit = tmp_path / f"explicit-{workflow_type}.duckdb"
    automatic = tmp_path / f"automatic-{workflow_type}.duckdb"
    assert VaultBuilder(renamed, explicit, workflow_type=workflow_type).build()[
        "validation"]["valid"]
    assert VaultBuilder(renamed / "results", automatic).build()[
        "validation"]["valid"]
    assert snapshot(original) == snapshot(explicit) == snapshot(automatic)
    with duckdb.connect(str(explicit), read_only=True) as connection:
        assert connection.execute(
            "SELECT value FROM metadata WHERE key='workflow_type'"
        ).fetchone()[0] == workflow_type
        assert connection.execute(
            "SELECT pipeline,pipeline_version FROM analysis"
        ).fetchone() == ("results", None)

    qc_bytes = qc.read_bytes()
    qc.unlink()
    layered = tmp_path / f"layered-{workflow_type}.duckdb"
    VaultBuilder(renamed, layered, workflow_type=workflow_type).build()
    qc.write_bytes(qc_bytes)
    assert VaultLayerAdder(layered, ["qc"]).add()["imported_files"] == 1
    assert snapshot(layered) == snapshot(original)


def test_renamed_layout_requires_explicit_type_when_ambiguous(pipeline, tmp_path):
    version = next(pipeline.iterdir())
    version.rename(pipeline / "results")
    workflow = next(pipeline.glob("*/*/main_workflow"))
    (workflow / "quantitate-ProteomicProtein-1.0.0").mkdir()
    with pytest.raises(ValueError, match="pass --workflow-type DIA or DDA"):
        VaultBuilder(pipeline, tmp_path / "ambiguous.duckdb").plan()
    assert VaultBuilder(pipeline, tmp_path / "dia.duckdb",
                        workflow_type="DIA").plan()["workflow_type"] == "DIA"


def test_renamed_dda_tmt_keeps_reporter_handling(tmp_path):
    root = dda_output(tmp_path, "TMT")
    next(root.iterdir()).rename(root / "results")
    channels = tmp_path / "channel.csv"
    channels.write_text("plex_id,channel,description\n"
                        "plex_01,126,reference\nplex_01,127N,treated\n")
    database = tmp_path / "tmt.duckdb"
    assert VaultBuilder(root, database, workflow_type="DDA",
                        tmt_channel_csv=channels).build()["validation"]["valid"]
    with duckdb.connect(str(database), read_only=True) as connection:
        assert connection.execute(
            "SELECT sample_name FROM samples ORDER BY sample_name"
        ).fetchall() == [("plex_01|126",), ("plex_01|127N",)]
        assert connection.execute("SELECT count(*) FROM tmt_runs").fetchone()[0] == 2


def test_explicit_workflow_type_checks_expected_artifacts(pipeline, tmp_path):
    with pytest.raises(ValueError, match="Missing DDA acquisition manifest"):
        VaultBuilder(pipeline, tmp_path / "wrong-dda.duckdb",
                     workflow_type="DDA").plan()
    dda = dda_output(tmp_path, "LFQ")
    with pytest.raises(ValueError, match="Expected one published DIA-NN"):
        VaultBuilder(dda, tmp_path / "wrong-dia.duckdb",
                     workflow_type="DIA").plan()


def test_existing_vault_without_workflow_type_can_add_layer(pipeline, tmp_path):
    database = tmp_path / "existing.duckdb"
    VaultBuilder(pipeline, database).build()
    with duckdb.connect(str(database)) as connection:
        connection.execute("DELETE FROM metadata WHERE key='workflow_type'")
    assert VaultLayerAdder(database, ["raw-search"]).add()["imported_files"] == 1
