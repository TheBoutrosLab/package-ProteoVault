"""Pure parsing and formatting operations; no file fixtures."""

import pytest

from proteovault.cli.inspect import VaultInspector
from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class TestTableFormatting:
    def test_format_table(self):
        table = VaultInspector.format_table(
            ("Name", "Count"), [("psms", "1"), ("protein_quant", "10")]
        )
        assert "Name" in table
        assert "protein_quant  10" in table
