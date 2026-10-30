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


class TestSourceImporters:
    def test_standalone_fasta_and_mzml_importers(self, tmp_path):
        from proteovault.database import ProteoVaultDatabase
        from proteovault.importers.fasta import FASTAImporter
        from proteovault.importers.mzml import MzMLImporter
        from proteovault.manifest import SourceManifest

        fasta = tmp_path / "proteins.fasta"
        shutil.copyfile(TEST_FIXTURES / "formats/proteins.fasta", fasta)
        mzml = tmp_path / "S1.mzML"
        shutil.copyfile(TEST_FIXTURES / "formats/S1.mzML", mzml)
        with ProteoVaultDatabase(tmp_path / "formats.duckdb") as database:
            database.initialize()
            manifest = SourceManifest(tmp_path)
            for path, name, cls in [
                (fasta, "fasta", FASTAImporter),
                (mzml, "mzml", MzMLImporter),
            ]:
                source = manifest.add_source(
                    database, dict(path=path, format=name, stage="input")
                )
                cls(path, database, source, "input").import_data()
            assert (
                database.connection.execute(
                    "SELECT rt_seconds FROM spectra"
                ).fetchone()[0]
                == 120
            )
            assert (
                database.connection.execute(
                    "SELECT sequence FROM proteins WHERE accession='P1'"
                ).fetchone()[0]
                == "PEPTIDE"
            )
            assert (
                database.connection.execute(
                    "SELECT count(*) FROM proteins WHERE is_decoy"
                ).fetchone()[0]
                == 1
            )
