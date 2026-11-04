"""Repair provenance after moving a completed analysis."""

import json
import shutil

from proteovault import VaultBuilder, VaultPathUpdater, VaultValidator
from proteovault.cli.main import main
from proteovault.database import ProteoVaultDatabase


def test_update_paths_after_move(pipeline, tmp_path, capsys):
    database = tmp_path / "results.duckdb"
    VaultBuilder(pipeline, database).build()
    old_root = pipeline
    new_root = tmp_path / "moved" / pipeline.name
    new_root.parent.mkdir()
    shutil.move(str(pipeline), str(new_root))

    args = ["update-paths", str(database), "--pipeline-output", str(new_root)]
    before = database.read_bytes()
    assert main(args + ["--dry-run"]) == 0
    assert json.loads(capsys.readouterr().out)["source_changes"]
    assert database.read_bytes() == before

    assert main(args) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["updated_sources"] > 0
    assert result["old_root"] == str(old_root)
    assert result["new_root"] == str(new_root)
    assert VaultValidator(database).validate()["valid"]
    with ProteoVaultDatabase(database, read_only=True) as vault:
        assert vault.connection.execute("SELECT root_path FROM analysis").fetchone()[0] == str(new_root)
        paths = vault.connection.execute(
            "SELECT path,relative_path FROM source_files WHERE format!='fasta'"
        ).fetchall()
        assert all(path == str(new_root / relative) for path, relative in paths)
    assert VaultPathUpdater(database, new_root).plan()["source_changes"] == []
