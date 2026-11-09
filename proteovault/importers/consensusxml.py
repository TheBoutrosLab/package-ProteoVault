"""Read protein inference groups from OpenMS consensusXML."""

import xml.etree.ElementTree as ET

from proteovault.importers.idxml import IdXMLImporter


class ConsensusXMLImporter(IdXMLImporter):
    """Import ProteinInferencer groups without treating consensus features as PSMs."""

    def import_data(self) -> dict[str, int]:
        self.counts = {"inference_groups": 0, "protein_hits": 0}
        refs = {}
        engine = "ProteinInferencer"
        self.identification_run_index = 0
        for event, element in ET.iterparse(self.path, events=("start", "end")):
            tag = element.tag.rsplit("}", 1)[-1]
            if event == "start":
                if tag == "IdentificationRun":
                    self.identification_run_index += 1
                    refs = {}
                continue
            if tag == "ProteinHit":
                accession = element.get("accession")
                if not accession:
                    raise ValueError(f"Missing protein accession: {self.path}")
                self.database.protein(accession, element.get("sequence") or None)
                refs[element.get("id")] = (None, accession)
                self.counts["protein_hits"] += 1
                element.clear()
            elif tag == "ProteinIdentification":
                parameters = [
                    (child.get("name", ""), child.get("value", ""))
                    for child in element
                    if child.tag.rsplit("}", 1)[-1] == "UserParam"
                    and child.get("name", "").startswith(
                        ("protein_group_", "indistinguishable_proteins_")
                    )
                ]
                self.import_groups(parameters, refs, engine, element.get("score_type"))
                element.clear()
        if not self.counts["protein_hits"]:
            raise ValueError(f"No protein identifications found: {self.path}")
        return self.record_counts()
