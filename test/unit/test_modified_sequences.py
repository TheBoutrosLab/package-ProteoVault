"""Pure parsing and formatting operations; no file fixtures."""

import pytest

from proteovault.cli.inspect import VaultInspector
from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class TestModifiedSequences:
    def test_unmodified_sequence(self):
        assert (
            IdXMLImporter.unmodified("M(Oxidation)PEPC(Carbamidomethyl)K") == "MPEPCK"
        )

    def test_unmodified_rejects_unbalanced_annotation(self):
        with pytest.raises(ValueError, match="Unbalanced"):
            IdXMLImporter.unmodified("M(Oxidation")
