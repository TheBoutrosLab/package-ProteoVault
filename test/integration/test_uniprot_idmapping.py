"""UniProt cross-references retain source provenance and original accessions."""

import gzip

import duckdb
import pytest

from proteovault import VaultBuilder, VaultInspector, VaultUniProtMappingAdder
from proteovault.cli.main import main


def prepare_pipeline(pipeline):
    for path in pipeline.rglob("*"):
        if path.suffix in {".idXML", ".tsv"}:
            content = path.read_text()
            if "P1" in content:
                path.write_text(content.replace("P1", "sp|O75083|WDR1_HUMAN"))


def mapping_file(tmp_path, content=None):
    path = tmp_path / "human.idmapping.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(content or (
            "O75083\tEnsembl\tENSG00000071127.19\n"
            "O75083\tEnsembl\tENSG00000071127.20\n"
            "O75083\tGene_Name\tWDR1\n"
            "O75083\tGene_Name\tWDR1\n"
            "O75083\tGeneID\t9948\n"
            "NOT_IN_VAULT\tGene_Name\tOTHER\n"
            "O75083\tEnsembl\t-\n"
        ))
    return path


def assert_mappings(database):
    with duckdb.connect(str(database), read_only=True) as connection:
        assert connection.execute(
            "SELECT p.accession, i.mapping_type, i.mapped_id "
            "FROM protein_identifiers i JOIN proteins p USING(protein_id) "
            "ORDER BY i.mapping_type, i.mapped_id"
        ).fetchall() == [
            ("sp|O75083|WDR1_HUMAN", "Ensembl", "ENSG00000071127.19"),
            ("sp|O75083|WDR1_HUMAN", "Ensembl", "ENSG00000071127.20"),
            ("sp|O75083|WDR1_HUMAN", "GeneID", "9948"),
            ("sp|O75083|WDR1_HUMAN", "Gene_Name", "WDR1"),
        ]
        assert connection.execute(
            "SELECT format, stage FROM source_files WHERE format='uniprot_idmapping'"
        ).fetchall() == [("uniprot_idmapping", "annotation")]
        assert connection.execute(
            "SELECT row_count FROM import_counts WHERE entity='protein_identifiers'"
        ).fetchall() == [(4,)]


def test_build_with_uniprot_mapping(pipeline, tmp_path):
    prepare_pipeline(pipeline)
    mapping = mapping_file(tmp_path)
    database = tmp_path / "mapped.duckdb"
    result = VaultBuilder(pipeline, database, uniprot_idmapping=mapping).build()
    assert result["validation"]["valid"]
    assert result["table_counts"]["protein_identifiers"] == 4
    assert "uniprot-idmapping" in VaultInspector(database).summary()["layers"]
    assert_mappings(database)


def test_add_mapping_is_atomic_and_idempotent(pipeline, tmp_path):
    prepare_pipeline(pipeline)
    database = tmp_path / "existing.duckdb"
    VaultBuilder(pipeline, database).build()
    mapping = mapping_file(tmp_path)
    assert main(["add-uniprot-idmapping", str(database), "--mapping", str(mapping), "--dry-run"]) == 0
    result = VaultUniProtMappingAdder(database, mapping).add()
    assert result["identifiers"] == 4
    assert result["matched_proteins"] == 1
    assert result["validation"]["valid"]
    assert_mappings(database)
    before = database.read_bytes()
    assert VaultUniProtMappingAdder(database, mapping).add()["status"] == "already_imported"
    assert database.read_bytes() == before
    bad = tmp_path / "bad.idmapping.gz"
    with gzip.open(bad, "wt") as handle:
        handle.write("O75083\tEnsembl\n")
    with pytest.raises(ValueError, match="expected three tab-separated columns"):
        VaultUniProtMappingAdder(database, bad).add()
    assert database.read_bytes() == before


def test_mapping_rejects_non_uniprot_proteins(pipeline, tmp_path):
    database = tmp_path / "unmapped.duckdb"
    mapping = mapping_file(tmp_path)
    with pytest.raises(ValueError, match="no unambiguous UniProt"):
        VaultBuilder(pipeline, database, uniprot_idmapping=mapping).build()
    assert not database.exists()
