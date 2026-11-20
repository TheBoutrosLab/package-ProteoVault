"""Select DIA or DDA from the workflow rather than its directory name."""

import json
import shutil
from pathlib import Path

import duckdb
import pytest

from proteovault import VaultBuilder, VaultLayerAdder
from proteovault.cli.main import main


def dda_output(tmp_path, mode):
    """Create the published DDA artifacts needed to exercise workflow routing."""
    root = tmp_path / "metapipeline-ProteomicDDA"
    workflow = root / "metapipeline-ProteomicDDA-1.0.0/TEST/main_workflow"
    manifest = workflow / "manifests/qc-mzml.csv"
    manifest.parent.mkdir(parents=True)
    if mode == "TMT":
        manifest.write_text(
            "sample_name,fraction_id,path\n"
            "plex_01,1,/data/fraction_1.mzML\n"
            "plex_01,2,/data/fraction_2.mzML\n"
        )
        (workflow / "quantitate-ProteomicTMT-1.0.0").mkdir()
        comments = (
            "# Files/samples associated with abundance values below: "
            "1: 'plex_01' ('tmt10plex_126'), "
            "2: 'plex_01' ('tmt10plex_127N')\n"
        )
    else:
        manifest.write_text(
            "sample_name,fraction_id,path\n"
            "S1,1,/data/S1.mzML\nS1,2,/data/S1_f2.mzML\n"
            "S2,1,/data/S2.mzML\n"
        )
        comments = (
            "# Files/samples associated with abundance values below: "
            "1: 'S1', 2: 'S2'\n"
        )
    output = (workflow / "quantitate-ProteomicProtein-1.0.0/TEST/"
              "OpenMS-3.5.0/output")
    output.mkdir(parents=True)
    (output / "TEST_ProteinQuantifier_protein.csv").write_text(
        comments + '"protein"\t"n_proteins"\t"protein_score"\t'
        '"n_peptides"\t"abundance_1"\t"abundance_2"\n'
        '"P1/P2"\t2\t0.9\t1\t100\t200\n'
    )
    (output / "TEST_ProteinQuantifier_peptide.csv").write_text(
        comments + '"peptide"\t"protein"\t"n_proteins"\t"charge"\t'
        '"abundance_1"\t"abundance_2"\t"fraction"\n'
        '"PEPTIDE(TMTpro)"\t"P1/P2"\t2\t0\t10\t20\t"all"\n'
    )
    source = (Path(__file__).parents[1] / "fixtures" /
              "metapipeline-ProteomicDIA-UniProtRP_ProteinProphet" /
              "metapipeline-ProteomicDIA-1.3.0/TEST/main_workflow" /
              "validate-ProteomicPSM-2.2.0/TEST/ProteinProphet/output" /
              "TEST-ProteinProphet.idXML")
    shutil.copyfile(source, output / "TEST_ProteinProphet_filtered.idXML")
    qc = workflow / "qc-ProteomicDDA-1.0.0/TEST/output/TEST_QC_metrics.tsv"
    qc.parent.mkdir(parents=True)
    qc_sample = "plex_01" if mode == "TMT" else "S1"
    qc.write_text(
        "dataset_id\tlayer\tsource_process\tsample_name\tsearch_engine\t"
        "fraction_id\tmetric_name\tgroup\tvalue\tunit\tstatus\tsource_file\n"
        f"TEST\tacquisition\tsummarize_mzML\t{qc_sample}\t\t1\t"
        "ms2_count\tall\t100\tcount\tavailable\tfraction_1.mzML\n"
        f"TEST\tacquisition\tsummarize_mzML\t{qc_sample}\t\t2\t"
        "ms2_count\treference\t120\tcount\tavailable\tfraction_2.mzML\n"
    )
    return root


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
