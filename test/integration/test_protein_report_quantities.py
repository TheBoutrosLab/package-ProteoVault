"""DIA-NN report protein quantities are source-specific group/run measurements."""

import csv

import duckdb
import pytest

from proteovault import VaultBuilder


def write_report(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "Protein.Group", "Stripped.Sequence", "Modified.Sequence",
            "Precursor.Charge", "Run", "Precursor.Quantity", "PG.Quantity",
            "PG.Normalised", "PG.MaxLFQ", "PG.TopN",
        ])
        writer.writerows(rows)


def test_report_protein_methods_are_distinct_and_deduplicated(pipeline, tmp_path):
    report = next(pipeline.rglob("*-report.tsv"))
    write_report(report, [
        ["P1;P2", "PEPTIDE", "PEPTIDE", 2, "S1", 19, 100, 110, 120, ""],
        ["P1;P2", "OTHER", "OTHER", 3, "S1", 20, 100, 110, 120, 130],
        ["P1;P2", "PEPTIDE", "PEPTIDE", 2, "S2", 21, "", 0, 300, ""],
    ])
    database = tmp_path / "quantities.duckdb"
    VaultBuilder(pipeline, database).build()

    with duckdb.connect(str(database), read_only=True) as connection:
        rows = connection.execute(
            "SELECT q.sample_name, q.quantification_method, q.quantity "
            "FROM protein_quant_long q JOIN source_files f USING(source_id) "
            "WHERE f.format='diann_report' "
            "ORDER BY q.sample_name, q.quantification_method"
        ).fetchall()
        assert rows == [
            ("S1", "PG.MaxLFQ", 120),
            ("S1", "PG.Normalised", 110),
            ("S1", "PG.Quantity", 100),
            ("S1", "PG.TopN", 130),
            ("S2", "PG.MaxLFQ", 300),
            ("S2", "PG.Normalised", 0),
            ("S2", "PG.Quantity", None),
            ("S2", "PG.TopN", None),
        ]
        assert connection.execute(
            "SELECT count(*) FROM protein_quant_readable WHERE method='matrix'"
        ).fetchone()[0] == 4
        assert connection.execute(
            "SELECT row_count FROM import_counts i JOIN source_files f USING(source_id) "
            "WHERE f.format='diann_report' AND i.entity='protein_quant'"
        ).fetchone()[0] == 8
        assert connection.execute(
            "SELECT count(*) FROM precursor_quant_readable "
            "WHERE method='Precursor.Quantity'"
        ).fetchone()[0] == 3
        assert connection.execute(
            "SELECT count(*) FROM protein_quant_readable WHERE method='PG.MaxLFQ'"
        ).fetchone()[0] == 2


def test_conflicting_repeated_report_quantity_rejects_core_build(pipeline, tmp_path):
    report = next(pipeline.rglob("*-report.tsv"))
    write_report(report, [
        ["P1;P2", "PEPTIDE", "PEPTIDE", 2, "S1", 19, 100, 110, 120, 130],
        ["P1;P2", "OTHER", "OTHER", 3, "S1", 20, 101, 110, 120, 130],
    ])
    database = tmp_path / "conflict.duckdb"
    with pytest.raises(ValueError, match="conflicting PG.Quantity"):
        VaultBuilder(pipeline, database).build()
    assert not database.exists()


def test_core_imports_all_available_precursor_report_methods_and_confidence(pipeline, tmp_path):
    report = next(pipeline.rglob("*-report.tsv"))
    report.write_text(
        report.read_text().replace("Precursor.Normalised\tQ.Value", "Precursor.Normalised\tPrecursor.Translated\tQ.Value")
        .replace("\tS1\t19\t20\t0.005", "\tS1\t19\t20\t21\t0.005")
    )
    database = tmp_path / "precursor_quantities.duckdb"
    VaultBuilder(pipeline, database).build()

    with duckdb.connect(str(database), read_only=True) as connection:
        assert connection.execute(
            "SELECT method, quantity, q_value, pep FROM precursor_quant_readable "
            "WHERE method != 'matrix' ORDER BY method"
        ).fetchall() == [
            ("Precursor.Normalised", 20, 0.005, 0.002),
            ("Precursor.Quantity", 19, 0.005, 0.002),
            ("Precursor.Translated", 21, 0.005, 0.002),
        ]
