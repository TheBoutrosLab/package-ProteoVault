"""Register DDA acquisition runs and TMT reporter channels."""

import csv
from pathlib import Path

from proteovault.importers.base import BaseImporter


TMT_REPORTERS = {
    "tmt6plex": ("126", "127", "128", "129", "130", "131"),
    "tmt10plex": ("126", "127N", "127C", "128N", "128C", "129N", "129C",
                  "130N", "130C", "131"),
    "tmt11plex": ("126", "127N", "127C", "128N", "128C", "129N", "129C",
                  "130N", "130C", "131N", "131C"),
    "tmt16plex": ("126", "127N", "127C", "128N", "128C", "129N", "129C",
                  "130N", "130C", "131N", "131C", "132N", "132C", "133N",
                  "133C", "134N"),
    "tmt18plex": ("126", "127N", "127C", "128N", "128C", "129N", "129C",
                  "130N", "130C", "131N", "131C", "132N", "132C", "133N",
                  "133C", "134N", "134C", "135N"),
}


class TMTChannelManifestImporter(BaseImporter):
    """Read the source channel CSV and identify one reference channel per plex."""

    def import_data(self) -> dict[str, int]:
        self.counts = {"tmt_plexes": 0, "tmt_channels": 0}
        references = {}
        declarations = {}
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != ["plex_id", "channel", "description"]:
                raise ValueError(f"Expected plex_id,channel,description: {self.path}")
            for row in reader:
                plex_name = row["plex_id"].strip()
                reporter = row["channel"].strip()
                description = row["description"].strip()
                if not plex_name or not reporter or "|" in plex_name + reporter:
                    raise ValueError(f"Invalid TMT plex or reporter label: {row}")
                key = (plex_name, reporter)
                if key in declarations:
                    raise ValueError(f"Duplicate TMT channel: {key}")
                declarations[key] = description or None
                if description.lower() == "reference":
                    if plex_name in references:
                        raise ValueError(f"Multiple reference channels for {plex_name}")
                    references[plex_name] = reporter
        for plex_name in sorted({plex for plex, _ in declarations}):
            supplied = [reporter for plex, reporter in declarations if plex == plex_name]
            reporters = TMT_REPORTERS.get(self.database.dda_plex_type, supplied)
            unknown = set(supplied) - set(reporters)
            if unknown:
                raise ValueError(
                    f"Reporter labels absent from {self.database.dda_plex_type}: "
                    f"{sorted(unknown)}"
                )
            plex = self.database.identity("tmt_plex", plex_name)
            self.database.dda_reporters[plex_name] = list(reporters)
            for reporter in reporters:
                channel = self.database.identity("tmt_channel", plex_name, reporter)
                sample = self.database.sample(f"{plex_name}|{reporter}")
                self.database.add("quant_channels", channel_id=channel,
                                  sample_id=sample, run_id=None, label="TMT")
                self.database.add(
                    "tmt_channels", channel_id=channel, plex_id=plex,
                    reporter_label=reporter,
                    channel_description=declarations.get((plex_name, reporter)),
                )
                self.database.dda_channels[(plex_name, reporter)] = channel
                self.database.dda_channel_descriptions[(plex_name, reporter)] = (
                    declarations.get((plex_name, reporter))
                )
                self.counts["tmt_channels"] += 1
            self.database.add(
                "tmt_plexes", plex_id=plex, plex_name=plex_name,
                reference_channel_id=(self.database.dda_channels[
                    (plex_name, references[plex_name])
                ] if plex_name in references else None),
            )
            self.counts["tmt_plexes"] += 1
        return self.record_counts()


class DDARunManifestImporter(BaseImporter):
    """Map each mzML fraction to an LFQ sample or TMT plex."""

    def import_data(self) -> dict[str, int]:
        self.counts = {"runs": 0, "tmt_runs": 0}
        seen = set()
        with self.path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            if columns == ["sample_name", "fraction_id", "path"]:
                fraction_column, path_column = "fraction_id", "path"
            elif columns == ["sample_name", "fraction", "mzML"]:
                fraction_column, path_column = "fraction", "mzML"
            else:
                raise ValueError(f"Unexpected DDA run manifest columns: {self.path}")
            for row in reader:
                name = row["sample_name"].strip()
                fraction = row[fraction_column].strip()
                path = row[path_column].strip()
                if not name or not fraction or not path or (name, fraction) in seen:
                    raise ValueError(f"Invalid or duplicate DDA run: {row}")
                seen.add((name, fraction))
                run = self.database.identity("dda_run", name, fraction)
                self.database.add(
                    "runs", run_id=run, run_name=Path(path).stem,
                    file_name=path,
                    sample_id=(self.database.sample(name)
                               if self.database.dda_mode == "LFQ" else None),
                )
                self.counts["runs"] += 1
                if self.database.dda_mode == "TMT":
                    if name not in self.database.dda_reporters:
                        raise ValueError(f"Unknown TMT plex in run manifest: {name}")
                    self.database.add("tmt_runs", run_id=run,
                                      plex_id=self.database.identity("tmt_plex", name))
                    self.counts["tmt_runs"] += 1
                elif name not in self.database.dda_channels:
                    channel = self.database.identity("dda_lfq_channel", name)
                    self.database.add("quant_channels", channel_id=channel,
                                      sample_id=self.database.sample(name),
                                      run_id=None, label="LFQ")
                    self.database.dda_channels[name] = channel
        if not seen:
            raise ValueError(f"Empty DDA run manifest: {self.path}")
        return self.record_counts()
