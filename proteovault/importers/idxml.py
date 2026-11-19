"""Stream OpenMS idXML while retaining source-specific PSM stages and group types."""

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

from proteovault.importers.base import BaseImporter


class IdXMLImporter(BaseImporter):
    """Read idXML without requiring pyOpenMS or discarding lower-ranked hits."""

    CONFIDENCE_SCORES = {
        "MS:1001491",
        "MS:1001493",
        "Posterior Error Probability",
        "Posterior Error Probability_score",
        "q-value",
        "PEP",
        "spectrum_q",
    }

    def __init__(self, *args, include_features: bool = True, **kwargs):
        super().__init__(*args, **kwargs)
        self.include_features = include_features

    @staticmethod
    def parameters(element):
        return {
            child.get("name"): child.get("value", "")
            for child in element
            if child.tag.rsplit("}", 1)[-1] == "UserParam"
        }

    @staticmethod
    def unmodified(sequence):
        # OpenMS modification names may contain nested parentheses.
        result = []
        depth = 0
        for character in sequence:
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
            elif depth == 0 and character.isupper():
                result.append(character)
        if depth != 0:
            raise ValueError(f"Unbalanced modification annotation: {sequence}")
        return "".join(result)

    def import_data(self) -> dict[str, int]:
        self.counts = {"psms": 0, "inference_groups": 0, "protein_hits": 0}
        dda_validation = (
            getattr(self.database, "is_dda", False) and self.stage == "validation"
        )
        self.run_id = (
            self.database.resolve_run(self.sample_name)
            if self.sample_name and not dda_validation else None
        )
        self.dda_runs = []
        self.spectra_data = []
        if dda_validation:
            self.database.flush()
            self.dda_runs = self.database.connection.execute(
                "SELECT r.run_name, r.run_id, s.sample_name, p.plex_name "
                "FROM runs r LEFT JOIN samples s ON r.sample_id=s.sample_id "
                "LEFT JOIN tmt_runs tr ON r.run_id=tr.run_id "
                "LEFT JOIN tmt_plexes p ON tr.plex_id=p.plex_id"
            ).fetchall()
        refs = {}
        engine = "unknown"
        protein_score_type = None
        identification_index = 0
        identification_run_index = 0
        group_parameters = []
        stack = []
        for event, element in ET.iterparse(self.path, events=("start", "end")):
            tag = element.tag.rsplit("}", 1)[-1]
            if event == "start":
                stack.append(element)
                if len(stack) == 1 and tag != "IdXML":
                    raise ValueError(f"Not an idXML document: {self.path}")
                if tag == "IdentificationRun":
                    identification_run_index += 1
                    self.identification_run_index = identification_run_index
                    refs = {}
                    group_parameters = []
                    self.spectra_data = []
                    engine = element.get("search_engine", "unknown")
                continue
            if tag == "UserParam" and element.get("name") == "spectra_data":
                parent_tag = (
                    stack[-2].tag.rsplit("}", 1)[-1] if len(stack) > 1 else None
                )
                if parent_tag in ("ProteinIdentification", "IdentificationRun"):
                    value = element.get("value", "").strip()
                    if value.startswith("[") and value.endswith("]"):
                        value = value[1:-1]
                    self.spectra_data = [
                        part.strip()
                        for part in next(csv.reader([value], skipinitialspace=True))
                    ]
            if tag == "ProteinHit":
                accession = element.get("accession", "")
                if not accession:
                    raise ValueError(f"Missing protein accession: {self.path}")
                protein = self.database.protein(
                    accession, element.get("sequence") or None
                )
                refs[element.get("id")] = (protein, accession)
                self.counts["protein_hits"] += 1
                element.clear()
            elif tag == "ProteinIdentification":
                protein_score_type = element.get("score_type")
                for child in element:
                    if child.tag.rsplit("}", 1)[-1] == "UserParam" and child.get(
                        "name", ""
                    ).startswith(("protein_group_", "indistinguishable_proteins_")):
                        group_parameters.append(
                            (child.get("name", ""), child.get("value", ""))
                        )
                self.import_groups(group_parameters, refs, engine, protein_score_type)
                element.clear()
            elif tag == "PeptideIdentification":
                identification_index += 1
                # Cohort inference idXML can refer to many source runs. Do not fabricate a run.
                if self.sample_name is not None or self.dda_runs:
                    self.import_identification(
                        element, refs, engine, identification_index
                    )
                element.clear()
                if identification_index % 5000 == 0:
                    self.database.flush()
            stack.pop()
            if stack and tag in (
                "ProteinHit",
                "ProteinIdentification",
                "PeptideIdentification",
            ):
                stack[-1].remove(element)
        if self.counts["protein_hits"] == 0 and identification_index == 0:
            raise ValueError(f"No identifications found in {self.path}")
        return self.record_counts()

    def import_groups(self, params, refs, engine, score_type):
        # OpenMS encodes groups as UserParam values: probability,PH_0,PH_1.
        for name, value in params:
            parts = [p.strip() for p in value.strip("[]").split(",")]
            if len(parts) < 2:
                raise ValueError(f"Unsupported idXML group encoding: {name}={value}")
            members = []
            for reference in parts[1:]:
                if reference not in refs:
                    raise ValueError(
                        f"Unknown protein reference {reference} in {self.path}"
                    )
                members.append(refs[reference][1])
            self.database.group(
                self.source_id,
                f"run{self.identification_run_index}:{name}",
                members,
                engine,
                name.rsplit("_", 1)[0],
                self.number(parts[0]),
                score_type,
            )
            self.counts["inference_groups"] += 1

    def import_identification(self, element, refs, engine, index):
        metadata = self.parameters(element)
        native = element.get("spectrum_reference") or metadata.get("spectrum_reference")
        if not native:
            raise ValueError(
                f"Missing spectrum reference in {self.path}, identification {index}"
            )
        run = self.run_id
        if self.dda_runs:
            run = self.resolve_dda_run(metadata, index, native)
        spectrum = self.database.identity("spectrum", run, native)
        self.database.add(
            "spectra",
            spectrum_id=spectrum,
            run_id=run,
            native_id=native,
            rt_seconds=self.number(element.get("RT")),
            ms_level=2,
            precursor_mz=self.number(element.get("MZ")),
            source_id=self.source_id,
        )
        score_type = element.get("score_type", "score")
        higher = element.get("higher_score_better", "false") == "true"
        hits = [h for h in element if h.tag.rsplit("}", 1)[-1] == "PeptideHit"]
        for hit_index, hit in enumerate(hits):
            modified = hit.get("sequence", "")
            sequence = self.unmodified(modified)
            precursor = self.database.precursor(
                sequence, modified, int(hit.get("charge", "0")), "OpenMS"
            )
            psm = self.database.identity("psm", self.source_id, index, hit_index)
            params = self.parameters(hit)
            self.database.add(
                "psms",
                psm_id=psm,
                spectrum_id=spectrum,
                precursor_id=precursor,
                search_engine=engine,
                candidate_rank=int(hit.get("rank", hit_index + 1)),
                target_decoy=params.get("target_decoy"),
                source_id=self.source_id,
                identification_index=index,
                hit_index=hit_index,
                stage=self.stage,
            )
            main_score_name = (
                "Posterior Error Probability"
                if score_type
                in (
                    "MS:1001493",
                    "Posterior Error Probability_score",
                    "Posterior Error Probability",
                    "PEP",
                )
                else score_type
            )
            self.database.add(
                "psm_scores",
                psm_id=psm,
                name=main_score_name,
                value=float(hit.get("score", "0")),
                higher_better=higher,
            )
            self.database.add(
                "score_aliases",
                score_id=self.database.identity("score_types", main_score_name)
                & ((1 << 31) - 1),
                source_name=score_type,
            )
            score_values = {main_score_name: float(hit.get("score", "0"))}
            for name, value in {**metadata, **params}.items():
                if self.include_features:
                    self.database.add(
                        "psm_metadata", psm_id=psm, name=name, value=value
                    )
                if not self.include_features and name not in self.CONFIDENCE_SCORES:
                    continue
                try:
                    numeric = self.number(value)
                except ValueError:
                    continue
                if numeric is not None:
                    if name in self.CONFIDENCE_SCORES and name != score_type:
                        canonical = (
                            "Posterior Error Probability"
                            if name
                            in (
                                "MS:1001493",
                                "Posterior Error Probability_score",
                                "Posterior Error Probability",
                                "PEP",
                            )
                            else name
                        )
                        if (
                            canonical in score_values
                            and score_values[canonical] != numeric
                        ):
                            canonical = name
                        if canonical not in score_values:
                            score_values[canonical] = numeric
                            self.database.add(
                                "psm_scores",
                                psm_id=psm,
                                name=canonical,
                                value=numeric,
                                higher_better=False,
                            )
                        self.database.add(
                            "score_aliases",
                            score_id=self.database.identity("score_types", canonical)
                            & ((1 << 31) - 1),
                            source_name=name,
                        )
                    if self.include_features:
                        self.database.add(
                            "psm_features", psm_id=psm, name=name, value=numeric
                        )
            for reference in hit.get("protein_refs", "").split():
                if reference not in refs:
                    raise ValueError(f"Unknown protein reference: {reference}")
                self.database.add(
                    "psm_proteins", psm_id=psm, protein_id=refs[reference][0]
                )
                self.database.add(
                    "peptide_proteins",
                    peptide_id=self.database.identity("peptide", sequence),
                    protein_id=refs[reference][0],
                )
            self.counts["psms"] += 1

    def resolve_dda_run(self, metadata, index, native):
        """Resolve a validated identification to one manifest mzML within its sample."""
        context = f"{self.path}, identification {index} ({native})"
        base = metadata.get("base_name", "").strip()
        if base:
            source = base
        else:
            if not self.spectra_data or any(not path for path in self.spectra_data):
                raise ValueError(f"Missing DDA spectra_data in {context}")
            merge_index = metadata.get("id_merge_index", "").strip()
            if len(self.spectra_data) == 1 and not merge_index:
                position = 0
            else:
                try:
                    position = int(merge_index)
                except ValueError as error:
                    raise ValueError(
                        f"Invalid DDA id_merge_index {merge_index!r} in {context}"
                    ) from error
                if position < 0 or position >= len(self.spectra_data):
                    raise ValueError(
                        f"Invalid DDA id_merge_index {position} in {context}"
                    )
            source = self.spectra_data[position]
        name = self.run_name(source)
        matches = [
            run_id for run_name, run_id, sample, plex in self.dda_runs
            if run_name == name and self.sample_name in (sample, plex)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one DDA run for {source!r} in sample/plex "
                f"{self.sample_name!r}; found {len(matches)} in {context}"
            )
        return matches[0]
