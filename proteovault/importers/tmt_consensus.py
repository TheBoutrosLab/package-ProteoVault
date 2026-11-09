"""Import identified TMT reporter observations from final consensusXML."""

from pathlib import Path
import xml.etree.ElementTree as ET

from proteovault.importers.base import BaseImporter
from proteovault.importers.idxml import IdXMLImporter


class TMTConsensusImporter(BaseImporter):
    """Retain one reporter value per source element and map."""

    def import_data(self) -> dict[str, int]:
        self.database.flush()
        format_name = self.database.connection.execute(
            "SELECT format FROM source_files WHERE source_id=?", [self.source_id]
        ).fetchone()[0]
        if format_name not in ("tmt_consensus_scan", "tmt_consensus_feature"):
            raise ValueError(f"Unsupported TMT consensus source format: {format_name}")
        scan_based = format_name == "tmt_consensus_scan"
        method = "TMT.ScanBased" if scan_based else "TMT.FeatureBased"
        self.counts = {"precursor_quant": 0}
        accepted = {}
        for run_name, native, precursor, spectrum in self.database.connection.execute(
            "SELECT r.run_name,s.native_id,p.precursor_id,s.spectrum_id "
            "FROM psms p JOIN spectra s USING(spectrum_id) JOIN runs r USING(run_id) "
            "JOIN processing_stages t USING(stage_id) "
            "WHERE t.name='validation' AND p.candidate_rank=1"
        ).fetchall():
            key = (run_name, native, precursor)
            if key in accepted and accepted[key] != spectrum:
                raise ValueError(f"Ambiguous accepted PSM for {key}")
            accepted[key] = spectrum
        maps = {}
        seen_observations = set()
        matched = 0
        for _, element in ET.iterparse(self.path, events=("end",)):
            tag = element.tag.rsplit("}", 1)[-1]
            if tag == "mapList":
                for entry in element:
                    if entry.tag.rsplit("}", 1)[-1] != "map":
                        continue
                    plex = entry.get("name")
                    label = entry.get("label", "")
                    reporter = label.rsplit("_", 1)[-1]
                    key = (plex, reporter)
                    if key not in self.database.dda_channels:
                        raise ValueError(f"Unknown TMT reporter map {key} in {self.path}")
                    maps[int(entry.get("id"))] = (plex, self.database.dda_channels[key])
                element.clear()
            elif tag == "consensusElement":
                identifications = element.findall("PeptideIdentification")
                if len(identifications) != 1:
                    element.clear()
                    continue
                identification = identifications[0]
                hits = identification.findall("PeptideHit")
                if len(hits) != 1:
                    element.clear()
                    continue
                hit = hits[0]
                modified = hit.get("sequence", "")
                charge = int(hit.get("charge", "0"))
                if not modified or charge <= 0:
                    raise ValueError(f"Invalid identified TMT feature {element.get('id')}")
                precursor = self.database.precursor(
                    IdXMLImporter.unmodified(modified), modified, charge, "OpenMS"
                )
                spectrum = None
                if scan_based:
                    metadata = IdXMLImporter.parameters(identification)
                    run_name = Path(metadata.get("base_name", "").replace("\\", "/")).name
                    native = identification.get("spectrum_reference")
                    scan_id = IdXMLImporter.parameters(element).get("scan_id")
                    if not native or native != scan_id:
                        raise ValueError(f"Scan identity mismatch in {self.path}: {element.get('id')}")
                    spectrum = accepted.get((run_name, native, precursor))
                    if spectrum is None:
                        element.clear()
                        continue
                    matched += 1
                centroid = element.find("centroid")
                rt = self.number(centroid.get("rt")) if centroid is not None else None
                source_element = element.get("id")
                if not source_element:
                    raise ValueError(f"Missing consensus element ID: {self.path}")
                for reporter in element.findall("./groupedElementList/element"):
                    map_id = int(reporter.get("map"))
                    if map_id not in maps:
                        raise ValueError(f"Unknown consensus map {map_id}: {self.path}")
                    observation = (source_element, map_id)
                    if observation in seen_observations:
                        raise ValueError(f"Repeated TMT source observation {observation}")
                    seen_observations.add(observation)
                    self.database.add(
                        "precursor_quant",
                        quant_id=self.database.identity("tmt_precursor_quant", self.source_id,
                                                        source_element, map_id),
                        precursor_id=precursor,
                        channel_id=maps[map_id][1],
                        source_id=self.source_id,
                        method=method,
                        quantity=self.number(reporter.get("it")),
                        source_element_id=source_element,
                        source_map_id=map_id,
                        rt_seconds=rt,
                        spectrum_id=spectrum,
                    )
                    self.counts["precursor_quant"] += 1
                element.clear()
                if self.counts["precursor_quant"] % 5000 == 0:
                    self.database.flush()
        if not maps:
            raise ValueError(f"No reporter maps in {self.path}")
        if scan_based and not matched:
            raise ValueError(f"No scan-based TMT elements matched accepted PSMs: {self.path}")
        return self.record_counts()
