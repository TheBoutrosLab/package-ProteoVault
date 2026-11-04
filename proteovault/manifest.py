"""Track imported artifacts and immutable source identities."""

from pathlib import Path

from proteovault.database import ProteoVaultDatabase


class SourceManifest:
    """Record source paths relative to the analysis root and import status."""

    def __init__(self, root: str | Path):
        self.root = Path(root).resolve()

    def add_source(self, database: ProteoVaultDatabase, artifact: dict) -> str:
        path = Path(artifact["path"])
        stat = path.stat()
        relative = (
            artifact.get("relative_path") or path.relative_to(self.root).as_posix()
        )
        identity = ("source", relative)
        if artifact["format"] == "diann_pg_report":
            identity += ("protein_quant",)
        key = database.identity(*identity)
        database.add(
            "source_files",
            source_id=key,
            path=str(path),
            relative_path=relative,
            format=artifact["format"],
            stage=artifact["stage"],
            size_bytes=stat.st_size,
            mtime_ns=stat.st_mtime_ns,
            status="imported",
            details=None,
        )
        return key
