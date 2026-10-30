"""Stream mzML spectrum metadata, excluding binary peak arrays."""

import xml.etree.ElementTree as ET

from proteovault.importers.base import BaseImporter


class MzMLImporter(BaseImporter):
    """Read scan metadata and convert scan start time to seconds."""

    def import_data(self):
        self.counts = {"spectra": 0}
        run = self.database.resolve_run(
            self.sample_name or self.run_name(self.path.name)
        )
        for event, element in ET.iterparse(self.path, events=("end",)):
            tag = element.tag.rsplit("}", 1)[-1]
            if tag == "binary":
                element.clear()
            if tag != "spectrum":
                continue
            native = element.get("id")
            if not native:
                raise ValueError(f"Missing spectrum ID: {self.path}")
            values = {}
            for child in element.iter():
                if child.tag.rsplit("}", 1)[-1] == "cvParam":
                    values[child.get("accession")] = child
            rt = values.get("MS:1000016")
            seconds = None
            if rt is not None:
                seconds = self.number(rt.get("value"))
                unit = rt.get("unitAccession")
                if unit == "UO:0000031":
                    seconds *= 60
                elif unit != "UO:0000010":
                    raise ValueError(f"Unsupported scan-time unit {unit}: {self.path}")

            def value(accession):
                item = values.get(accession)
                return self.number(item.get("value")) if item is not None else None

            self.database.add(
                "spectra",
                spectrum_id=self.database.identity("spectrum", run, native),
                run_id=run,
                native_id=native,
                rt_seconds=seconds,
                ms_level=value("MS:1000511"),
                precursor_mz=value("MS:1000744"),
                precursor_charge=value("MS:1000041"),
                source_id=self.source_id,
            )
            self.counts["spectra"] += 1
            element.clear()
            if self.counts["spectra"] % 2000 == 0:
                self.database.flush()
        if not self.counts["spectra"]:
            raise ValueError(f"No spectra found: {self.path}")
        return self.record_counts()
