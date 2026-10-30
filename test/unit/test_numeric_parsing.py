"""Pure parsing and formatting operations; no file fixtures."""

import pytest

from proteovault.cli.inspect import VaultInspector
from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class TestNumericParsing:
    @pytest.mark.parametrize(
        ("value", "expected"),
        [(None, None), ("NA", None), ("  ", None), ("0", 0.0), ("1.25", 1.25)],
    )
    def test_number(self, value, expected):
        assert BaseImporter.number(value) == expected

    def test_number_rejects_nonfinite(self):
        with pytest.raises(ValueError, match="Nonfinite"):
            BaseImporter.number("inf")
