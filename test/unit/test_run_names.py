"""Pure parsing and formatting operations; no file fixtures."""

import pytest

from proteovault.cli.inspect import VaultInspector
from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class TestRunNames:
    @pytest.mark.parametrize(
        ("path", "expected"),
        [("/data/S1.mzML", "S1"), (r"C:\\data\\S2.raw", "S2"), ("sample.d", "sample")],
    )
    def test_run_name(self, path, expected):
        assert BaseImporter.run_name(path) == expected
