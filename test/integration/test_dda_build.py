"""Build DDA LFQ and TMT vaults from representative published artifacts."""

import shutil
import json
from pathlib import Path

import duckdb
import pytest

from proteovault import VaultBuilder


INFERENCE = (Path(__file__).parents[1] / "fixtures" /
             "metapipeline-ProteomicDIA-UniProtRP_ProteinProphet" /
             "metapipeline-ProteomicDIA-1.3.0/TEST/main_workflow" /
             "validate-ProteomicPSM-2.2.0/TEST/ProteinProphet/output" /
             "TEST-ProteinProphet.idXML")


def dda_output(tmp_path, mode):
    root = tmp_path / "metapipeline-ProteomicDDA"
    workflow = root / "metapipeline-ProteomicDDA-1.0.0/TEST/main_workflow"
    manifest = workflow / "manifests/qc-mzml.csv"
    manifest.parent.mkdir(parents=True)
    if mode == "TMT":
        manifest.write_text(
            "sample_name,fraction_id,path\nplex_01,1,/data/fraction_1.mzML\n"
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
            "sample_name,fraction_id,path\nS1,1,/data/S1.mzML\n"
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
    shutil.copy(INFERENCE, output / "TEST_ProteinProphet_filtered.idXML")
    qc = workflow / "qc-ProteomicDDA-1.0.0/TEST/output/TEST_QC_metrics.tsv"
    qc.parent.mkdir(parents=True)
    qc.write_text(
        "dataset_id\tacquisition_mode\tlayer\tsample_name\tfraction_id\t"
        "search_engine\tinference_route\tsource_process\tmetric_name\tgroup\t"
        "value\tunit\tstatus\tsource_file\tcalculation_version\n"
        f"TEST\tDDA\tquantification\t{'plex_01' if mode == 'TMT' else 'S1'}\t\t"
        "\tProteinProphet\tsummarize_ProteinQuantifier\t"
        "quantified_feature_count\tall\t1\tcount\tavailable\tprotein.csv\t1.0.0\n"
    )
    return root


@pytest.mark.parametrize("mode", ["LFQ", "TMT"])
def test_dda_quantitation_and_channel_identity(tmp_path, mode):
    root = dda_output(tmp_path, mode)
    channels = tmp_path / "channel.csv"
    if mode == "TMT":
        channels.write_text(
            "plex_id,channel,description\n"
            "plex_01,126,reference\nplex_01,127N,treated\n"
        )
    output = tmp_path / "dda.duckdb"
    result = VaultBuilder(root, output,
                          tmt_channel_csv=channels if mode == "TMT" else None).build()
    assert result["validation"]["valid"]
    if mode == "TMT":
        assert any("has no OpenMS output" in warning
                   for warning in result["validation"]["warnings"])
    with duckdb.connect(str(output), read_only=True) as connection:
        names = [row[0] for row in connection.execute(
            "SELECT sample_name FROM samples ORDER BY sample_name"
        ).fetchall()]
        assert names == (["plex_01|126", "plex_01|127N"] if mode == "TMT"
                         else ["S1", "S2"])
        assert connection.execute("SELECT count(*) FROM protein_quant").fetchone()[0] == 2
        assert connection.execute("SELECT count(*) FROM peptide_quant").fetchone()[0] == 2
        assert connection.execute(
            "SELECT acquisition_mode, inference_route FROM qc_metrics"
        ).fetchone() == ("DDA", "ProteinProphet")
        assert connection.execute(
            "SELECT charge FROM peptide_quant_features"
        ).fetchone()[0] == 0
        assert connection.execute(
            "SELECT p.sequence, f.source_peptide FROM peptide_quant_features f "
            "JOIN peptides p USING(peptide_id)"
        ).fetchone() == ("PEPTIDE", "PEPTIDE(TMTpro)")
        if mode == "TMT":
            assert connection.execute("SELECT count(*) FROM tmt_runs").fetchone()[0] == 2
            assert connection.execute(
                "SELECT s.sample_name FROM tmt_plexes p "
                "JOIN quant_channels c ON p.reference_channel_id=c.channel_id "
                "JOIN samples s USING(sample_id)"
            ).fetchone()[0] == "plex_01|126"


def test_tmt_detailed_abundance_headers(tmp_path):
    root = dda_output(tmp_path, "TMT")
    output_dir = next(root.rglob("OpenMS-3.5.0/output"))
    for suffix in ("protein", "peptide"):
        path = output_dir / f"TEST_ProteinQuantifier_{suffix}.csv"
        text = path.read_text()
        text = text[text.index('"' + suffix + '"'):]
        text = text.replace('"abundance_1"', '"abundance|plex_01|ch1"')
        text = text.replace('"abundance_2"', '"abundance|plex_01|ch2"')
        path.write_text(text)
    channels = tmp_path / "channel.csv"
    channels.write_text(
        "plex_id,channel,description\n"
        "plex_01,126,reference\nplex_01,127N,treated\n"
    )
    output = tmp_path / "detailed.duckdb"
    assert VaultBuilder(root, output, tmt_channel_csv=channels).build()[
        "validation"
    ]["valid"]


def test_tmt_numbered_sample_headers_without_reference(tmp_path):
    root = dda_output(tmp_path, "TMT")
    workflow = next(root.rglob("main_workflow"))
    (workflow / "manifests/qc-mzml.csv").rename(
        workflow / "manifests/validation-mzml.csv"
    )
    output_dir = next(root.rglob("OpenMS-3.5.0/output"))
    for suffix in ("protein", "peptide"):
        path = output_dir / f"TEST_ProteinQuantifier_{suffix}.csv"
        content = path.read_text()
        content = content[content.index('"' + suffix + '"'):]
        content = content.replace('"abundance_1"', '"abundance_sample1[unknown]"')
        content = content.replace('"abundance_2"', '"abundance_sample2[unknown]"')
        path.write_text(content)
    channels = tmp_path / "channel.csv"
    channels.write_text(
        "plex_id,channel,description\n"
        "plex_01,126,control\nplex_01,127N,treated\n"
    )
    output = tmp_path / "numbered.duckdb"
    assert VaultBuilder(root, output, tmt_channel_csv=channels).build()[
        "validation"
    ]["valid"]
    with duckdb.connect(str(output), read_only=True) as connection:
        assert connection.execute(
            "SELECT reference_channel_id FROM tmt_plexes"
        ).fetchone()[0] is None
        assert connection.execute(
            "SELECT sample_name FROM samples ORDER BY sample_name"
        ).fetchall() == [("plex_01|126",), ("plex_01|127N",)]


def test_tmt_omitted_reporter_is_retained_without_shifting_columns(tmp_path):
    root = dda_output(tmp_path, "TMT")
    (root / "run.resolved.yaml").write_text(json.dumps({
        "quantitate_ProteomicTMT": {"plex_type": "tmt10plex"}
    }))
    output_dir = next(root.rglob("OpenMS-3.5.0/output"))
    for suffix in ("protein", "peptide"):
        path = output_dir / f"TEST_ProteinQuantifier_{suffix}.csv"
        rows = path.read_text().splitlines()
        header = rows[-2].split("\t")
        values = rows[-1].split("\t")
        header = header[:4] + [f'"abundance_sample{i}[unknown]"'
                               for i in range(1, 11)]
        if suffix == "protein":
            values = values[:4] + [str(i * 10) for i in range(1, 11)]
        else:
            header += ['"fraction"']
            values = values[:4] + [str(i * 10) for i in range(1, 11)] + ['"all"']
        path.write_text("\t".join(header) + "\n" + "\t".join(values) + "\n")
    channels = tmp_path / "channel.csv"
    reporters = ("126", "127N", "127C", "128N", "128C", "129N", "129C",
                 "130N", "130C")
    channels.write_text("plex_id,channel,description\n" + "".join(
        f"plex_01,{reporter},sample_{index}\n"
        for index, reporter in enumerate(reporters, start=1)
    ))
    output = tmp_path / "omitted.duckdb"
    assert VaultBuilder(root, output, tmt_channel_csv=channels).build()[
        "validation"
    ]["valid"]
    with duckdb.connect(str(output), read_only=True) as connection:
        assert connection.execute("SELECT count(*) FROM tmt_channels").fetchone()[0] == 10
        assert connection.execute(
            "SELECT channel_description FROM tmt_channels WHERE reporter_label='131'"
        ).fetchone()[0] is None
        assert connection.execute(
            "SELECT q.quantity FROM protein_quant_long q "
            "WHERE q.sample_name='plex_01|131'"
        ).fetchone()[0] == 100


@pytest.mark.parametrize("mode", ["ScanBased", "FeatureBased"])
def test_tmt_consensus_preserves_repeated_precursor_observations(tmp_path, mode):
    root = dda_output(tmp_path, "TMT")
    workflow = next(root.rglob("main_workflow"))
    (workflow / "quantitate-ProteomicTMT-1.0.0").rmdir()
    filtered = (workflow / "validate-ProteomicPSM-1.0.0/plex_01/"
                "OpenMS-3.5.0/output/plex_01_filter.idXML")
    filtered.parent.mkdir(parents=True)
    filtered.write_text('''<IdXML><IdentificationRun search_engine="MSFragger">
      <PeptideIdentification spectrum_reference="scan=1" RT="100" MZ="500"
        score_type="q-value"><UserParam name="base_name" value="fraction_1"/>
        <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2" score="0.001"/></PeptideIdentification>
      <PeptideIdentification spectrum_reference="scan=2" RT="200" MZ="500"
        score_type="q-value"><UserParam name="base_name" value="fraction_2"/>
        <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2" score="0.002"/></PeptideIdentification>
    </IdentificationRun></IdXML>''')
    output_dir = (workflow / "quantitate-ProteomicTMT-2.0.0/"
                  "OpenMS-3.5.0/output")
    output_dir.mkdir(parents=True)
    (output_dir / f"TEST_{mode}.consensusXML").write_text("<consensusXML/>")
    features = []
    for index, rt in enumerate((100, 200), start=1):
        spectrum = f'scan={index}'
        base_name = f'fraction_{index}'
        scan = (f'<UserParam name="scan_id" value="{spectrum}"/>'
                if mode == "ScanBased" else "")
        features.append(f'''<consensusElement id="e_{index}" charge="2">
          <centroid rt="{rt}" mz="500" it="1000"/>
          <groupedElementList>
            <element map="0" it="{index * 10}"/>
            <element map="1" it="{index * 20}"/>
          </groupedElementList>{scan}
          <PeptideIdentification spectrum_reference="{spectrum}">
            <UserParam name="base_name" value="{base_name}"/>
            <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2"/>
            <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2"/>
          </PeptideIdentification>
          <PeptideIdentification spectrum_reference="{spectrum}">
            <UserParam name="base_name" value="{base_name}"/>
            <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2"/>
          </PeptideIdentification>
        </consensusElement>''')
    features.extend([
        '<consensusElement id="no_id"><groupedElementList>'
        '<element map="0" it="10"/></groupedElementList></consensusElement>',
        '<consensusElement id="no_hit"><PeptideIdentification/>'
        '<groupedElementList><element map="0" it="10"/>'
        '</groupedElementList></consensusElement>',
    ])
    if mode == "ScanBased":
        features.append('''<consensusElement id="unmatched">
          <UserParam name="scan_id" value="scan=3"/>
          <PeptideIdentification spectrum_reference="scan=3">
            <UserParam name="base_name" value="fraction_1"/>
            <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2"/>
          </PeptideIdentification></consensusElement>''')
    else:
        features.append('''<consensusElement id="ambiguous">
          <PeptideIdentification>
            <PeptideHit sequence="PEPTIDE(TMTpro)" charge="2"/>
            <PeptideHit sequence="OTHER(TMTpro)" charge="2"/>
          </PeptideIdentification></consensusElement>''')
    (output_dir / "TEST_IDConflictResolver.consensusXML").write_text(
        '<consensusXML><mapList><map id="0" name="plex_01" label="tmt10plex_126"/>'
        '<map id="1" name="plex_01" label="tmt10plex_127N"/></mapList>'
        + "".join(features) + '</consensusXML>')
    channels = tmp_path / "channel.csv"
    channels.write_text("plex_id,channel,description\n"
                        "plex_01,126,reference\nplex_01,127N,treated\n")
    output = tmp_path / f"{mode}.duckdb"
    result = VaultBuilder(root, output, tmt_channel_csv=channels).build()
    assert result["validation"]["valid"]
    warnings = result["validation"]["warnings"]
    assert any("no_identification" in warning for warning in warnings)
    assert any("no_hit" in warning for warning in warnings)
    assert any(("no_accepted_psm" if mode == "ScanBased" else
                "ambiguous_identification") in warning for warning in warnings)
    with duckdb.connect(str(output), read_only=True) as connection:
        rows = connection.execute(
            "SELECT sample_name,quantity,rt_seconds,spectrum_id,source_element_id "
            "FROM precursor_quant_long ORDER BY sample_name,rt_seconds"
        ).fetchall()
        assert len(rows) == 4
        assert [row[1] for row in rows] == [10, 20, 20, 40]
        assert [row[2] for row in rows] == [100, 200, 100, 200]
        assert len({row[4] for row in rows}) == 2
        assert all(row[3] is not None for row in rows) == (mode == "ScanBased")


@pytest.mark.parametrize("explicit_protein_headers", [False, True])
def test_two_plex_protein_columns_require_explicit_identity(tmp_path, explicit_protein_headers):
    root = dda_output(tmp_path, "TMT")
    workflow = next(root.rglob("main_workflow"))
    manifest = workflow / "manifests/qc-mzml.csv"
    manifest.write_text(manifest.read_text().replace(
        "plex_01,2,/data/fraction_2.mzML", "plex_02,1,/data/fraction_2.mzML"))
    output_dir = next(root.rglob("OpenMS-3.5.0/output"))
    protein = output_dir / "TEST_ProteinQuantifier_protein.csv"
    peptide = output_dir / "TEST_ProteinQuantifier_peptide.csv"
    if explicit_protein_headers:
        protein.write_text(protein.read_text().replace('"abundance_1"',
            '"abundance|plex_01|ch1"').replace('"abundance_2"',
            '"abundance|plex_02|ch1"'))
    else:
        protein.write_text(protein.read_text().replace('"abundance_1"',
            '"abundance_sample1[unknown]"').replace('"abundance_2"',
            '"abundance_sample2[unknown]"'))
    peptide.write_text(peptide.read_text().replace('"abundance_1"',
        '"abundance|plex_01|ch1"').replace('"abundance_2"',
        '"abundance|plex_02|ch1"'))
    channels = tmp_path / "channel.csv"
    channels.write_text("plex_id,channel,description\n"
                        "plex_01,126,reference\nplex_02,126,reference\n")
    output = tmp_path / "two_plex.duckdb"
    if not explicit_protein_headers:
        with pytest.raises(ValueError, match="numbered protein columns do not identify"):
            VaultBuilder(root, output, tmt_channel_csv=channels).build()
        assert not output.exists()
        return
    assert VaultBuilder(root, output, tmt_channel_csv=channels).build()[
        "validation"]["valid"]
    with duckdb.connect(str(output), read_only=True) as connection:
        assert connection.execute(
            "SELECT sample_name,quantity FROM protein_quant_long ORDER BY sample_name"
        ).fetchall() == [("plex_01|126", 100), ("plex_02|126", 200)]
