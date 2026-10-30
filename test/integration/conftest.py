"""Prepared pipeline files copied for each integration test."""

import shutil
from pathlib import Path

import pytest


TEST_FIXTURES = Path(__file__).parents[1] / "fixtures"


@pytest.fixture
def pipeline(tmp_path):
    source = TEST_FIXTURES / "metapipeline-ProteomicDIA-UniProtRP_ProteinProphet"
    destination = tmp_path / source.name
    shutil.copytree(source, destination)
    return destination
