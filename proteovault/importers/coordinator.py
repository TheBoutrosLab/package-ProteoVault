"""Shared orchestration for building a vault and adding optional source layers."""

import logging
import json
from time import perf_counter

from proteovault.importers.diann import DIANNImporter
from proteovault.importers.fasta import FASTAImporter
from proteovault.importers.idxml import IdXMLImporter
from proteovault.importers.pin import PINImporter
from proteovault.importers.qc import QCImporter
from proteovault.importers.uniprot_idmapping import UniProtIdMappingImporter
from proteovault.importers.dda_manifests import (
    DDARunManifestImporter, TMTChannelManifestImporter, TMT_REPORTERS,
)
from proteovault.importers.openms_quant import OpenMSQuantImporter
from proteovault.importers.consensusxml import ConsensusXMLImporter
from proteovault.importers.tmt_consensus import TMTConsensusImporter
from proteovault.manifest import SourceManifest
from proteovault.discovery import PipelineOutputDiscovery
from proteovault.sample_manifest import SampleManifest

logger = logging.getLogger(__name__)


class ImportCoordinator:
    """Register sources and dispatch format importers within the caller's transaction."""

    def __init__(self, database, root, pipeline, acquisition_mode):
        self.database = database
        self.manifest = SourceManifest(root)
        workflow, _ = PipelineOutputDiscovery(root).analysis()
        self.database.is_dda = pipeline == "metapipeline-ProteomicDDA"
        self.database.dda_mode = acquisition_mode
        self.database.dda_channels = {}
        self.database.dda_reporters = {}
        self.database.dda_plex_type = None
        if self.database.dda_mode == "TMT":
            resolved = workflow.parents[2] / "run.resolved.yaml"
            if resolved.is_file():
                settings = json.loads(resolved.read_text(encoding="utf-8"))
                plex_type = settings.get("quantitate_ProteomicTMT", {}).get("plex_type")
                if plex_type not in TMT_REPORTERS:
                    raise ValueError(f"Unsupported TMT plex type: {plex_type}")
                self.database.dda_plex_type = plex_type
        manifest_path = workflow / "manifests" / "quantitation.csv"
        if pipeline == "metapipeline-ProteomicDIA" and manifest_path.is_file():
            self.database.sample_manifest = SampleManifest(manifest_path)

    def import_artifacts(self, artifacts):
        warnings = []
        if self.database.sample_manifest:
            for name in sorted(self.database.sample_manifest.samples):
                self.database.sample(name)
        for index, artifact in enumerate(artifacts, start=1):
            started = perf_counter()
            logger.info(
                "Importing file %d/%d [%s, %s]: %s",
                index,
                len(artifacts),
                artifact["stage"],
                artifact["format"],
                artifact["path"],
            )
            source = self.manifest.add_source(self.database, artifact)
            options = {}
            if artifact["format"] == "idxml":
                options["include_features"] = artifact["stage"] == "search"
            if artifact["format"] == "fasta":
                options["only_existing"] = True
            importer = self.importer_class(artifact["format"])(
                path=artifact["path"],
                database=self.database,
                source_id=source,
                stage=artifact["stage"],
                sample_name=artifact["sample_name"],
                **options,
            )
            counts = importer.import_data()
            warnings.extend(getattr(importer, "warnings", []))
            logger.info(
                "Imported %s in %.2f s: %s",
                artifact["path"],
                perf_counter() - started,
                ", ".join(f"{name}={count}" for name, count in counts.items()),
            )
        self.database.flush()
        return warnings

    @staticmethod
    def importer_class(format_name):
        if format_name.startswith("diann_"):
            return DIANNImporter
        return {
            "idxml": IdXMLImporter,
            "pin": PINImporter,
            "qc_metrics": QCImporter,
            "uniprot_idmapping": UniProtIdMappingImporter,
            "fasta": FASTAImporter,
            "dda_run_manifest": DDARunManifestImporter,
            "tmt_channel_manifest": TMTChannelManifestImporter,
            "openms_protein_quant": OpenMSQuantImporter,
            "openms_peptide_quant": OpenMSQuantImporter,
            "consensusxml": ConsensusXMLImporter,
            "tmt_consensus_scan": TMTConsensusImporter,
            "tmt_consensus_feature": TMTConsensusImporter,
        }[format_name]
