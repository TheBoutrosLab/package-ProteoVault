"""Discover the core or optional layers of one published proteomics analysis."""

import re
from pathlib import Path


class PipelineOutputDiscovery:
    """Resolve one analysis and select artifacts for an explicitly requested layer."""

    LAYERS = ("core", "raw-search", "msbooster", "qc")

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve()

    def analysis(self) -> tuple[Path, dict]:
        if not self.root.is_dir():
            raise FileNotFoundError(
                f"Pipeline output directory does not exist: {self.root}"
            )
        pattern = r"metapipeline-Proteomic(?:DIA|DDA)-[0-9].*"
        if re.fullmatch(pattern, self.root.name):
            candidates = [self.root]
        else:
            candidates = sorted(
                p for p in self.root.glob("metapipeline-ProteomicD*-*")
                if p.is_dir() and re.fullmatch(pattern, p.name)
            )
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one DDA or DIA metapipeline version under {self.root}; found {len(candidates)}"
            )
        version_root = candidates[0]
        workflows = sorted(version_root.glob("*/main_workflow"))
        if len(workflows) != 1:
            raise ValueError(
                f"Expected one dataset/main_workflow under {version_root}; found {len(workflows)}"
            )
        workflow = workflows[0]
        match = re.fullmatch(r"(metapipeline-ProteomicD(?:IA|DA))-(.+)", version_root.name)
        assert match is not None
        pipeline, version = match.groups()
        return workflow, dict(
            dataset_id=workflow.parent.name,
            pipeline=pipeline,
            pipeline_version=version,
            root=str(self.root),
        )

    def discover(self, layer: str = "core") -> dict:
        if layer not in self.LAYERS:
            raise ValueError(f"Unknown layer: {layer}")
        workflow, analysis = self.analysis()
        if analysis["pipeline"] == "metapipeline-ProteomicDDA":
            return self.discover_dda(workflow, analysis, layer)
        artifacts = []
        warnings = []
        if layer == "core":
            quant_roots = sorted(
                workflow.glob("quantitate-ProteomicDIA-*/*/output/DIA-NN")
            )
            if len(quant_roots) != 1:
                raise ValueError(
                    f"Expected one published DIA-NN output under {workflow}; found {len(quant_roots)}"
                )
            quant = quant_roots[0]
            if len(list(quant.glob("*-report.pg_matrix.tsv"))) != 1:
                raise ValueError(
                    f"Expected one DIA-NN protein-group matrix under {quant}"
                )
            for pattern, format_name in [
                ("*-report.pg_matrix.tsv", "diann_pg_matrix"),
                ("*-report.pr_matrix.tsv", "diann_pr_matrix"),
                ("*-report.tsv", "diann_report"),
            ]:
                matches = sorted(quant.glob(pattern))
                if len(matches) > 1:
                    raise ValueError(f"Ambiguous {format_name}: {matches}")
                artifacts.extend(
                    self.artifact(p, format_name, "quantitation") for p in matches
                )
                if format_name == "diann_report" and not matches:
                    warnings.append("No DIA-NN report quantities found")
            validation = self.child_pipeline(workflow, "validate-ProteomicPSM-*")
            if validation:
                artifacts.extend(
                    self.artifact(p, "idxml", "validation", p.parents[2].name)
                    for p in sorted(validation.glob("*/OpenMS-*/output/*_filter.idXML"))
                )
                artifacts.extend(
                    self.artifact(p, "idxml", "inference")
                    for p in sorted(validation.glob("*/*/output/*.idXML"))
                    if p.parents[1].name in ("ProteinProphet", "Fido", "Epifany")
                )
            artifacts.extend(self.qc_artifacts(workflow, analysis["dataset_id"]))
            if not any(item["stage"] == "validation" for item in artifacts):
                warnings.append(
                    "No supported FDR-filtered sample idXML found; validated PSMs will be empty"
                )
            if not any(item["stage"] == "inference" for item in artifacts):
                warnings.append(
                    "No supported cohort inference idXML found; quantitative groups remain available"
                )
            if not any(item["format"] in ("diann_pr_matrix",) for item in artifacts):
                warnings.append("No precursor quantities found")
        elif layer == "qc":
            artifacts.extend(self.qc_artifacts(workflow, analysis["dataset_id"]))
        elif layer == "raw-search":
            search = self.child_pipeline(workflow, "search-ProteomicDIA-*")
            if search:
                artifacts.extend(
                    self.artifact(
                        p, "idxml", "search", re.sub(r"_MSFragger$", "", p.stem)
                    )
                    for p in sorted(search.glob("*/output/MSFragger/*.idXML"))
                )
        elif layer == "msbooster":
            validation = self.child_pipeline(workflow, "validate-ProteomicPSM-*")
            if validation:
                artifacts.extend(
                    self.artifact(p, "pin", "MSBooster", p.parents[3].name)
                    for p in sorted(
                        validation.glob("*/MSBooster/output/engine_*/*_edited.pin")
                    )
                )
        if layer != "core" and not artifacts:
            raise ValueError(
                f"No supported files found for layer {layer} under {workflow}"
            )
        return {
            **analysis,
            "layer": layer,
            "artifacts": artifacts,
            "warnings": warnings,
        }

    def discover_dda(self, workflow: Path, analysis: dict, layer: str) -> dict:
        """Select published DDA results without treating them as DIA-NN files."""
        dataset = analysis["dataset_id"]
        artifacts = []
        warnings = []
        tmt_directories = sorted(workflow.glob("quantitate-ProteomicTMT-*"))
        mode = "TMT" if tmt_directories else "LFQ"
        if len(tmt_directories) > 1:
            raise ValueError(f"Expected one TMT quantitation directory under {workflow}")
        protein_roots = sorted(workflow.glob("quantitate-ProteomicProtein-*/" + dataset + "/OpenMS-*/output"))
        if layer == "core":
            run_manifest = workflow / "manifests" / "qc-mzml.csv"
            if not run_manifest.is_file():
                run_manifest = workflow / "manifests" / "validation-mzml.csv"
            if not run_manifest.is_file():
                raise ValueError(f"Missing DDA acquisition manifest under {workflow}")
            artifacts.append(self.artifact(run_manifest, "dda_run_manifest", "design"))
            validation = self.child_pipeline(workflow, "validate-ProteomicPSM-*")
            if validation:
                artifacts.extend(self.artifact(path, "idxml", "validation")
                                 for path in sorted(validation.glob(
                                     "*/OpenMS-*/output/*_filter.idXML")))
            if len(protein_roots) != 1:
                raise ValueError(
                    f"Expected one DDA protein quantitation output under {workflow}; found {len(protein_roots)}"
                )
            root = protein_roots[0]
            # Import peptide quantities before proteins; keep this order explicit.
            # Each file must map its own abundance headers to channels.
            for suffix, format_name in (
                ("_ProteinQuantifier_peptide.csv", "openms_peptide_quant"),
                ("_ProteinQuantifier_protein.csv", "openms_protein_quant"),
            ):
                matches = sorted(root.glob("*" + suffix))
                if len(matches) != 1:
                    raise ValueError(f"Expected one {format_name} under {root}; found {len(matches)}")
                artifacts.append(self.artifact(matches[0], format_name, "quantitation"))
            proteinprophet = sorted(root.glob("*_ProteinProphet_filtered.idXML"))
            if not proteinprophet:
                proteinprophet = sorted(
                    root.parent.glob("intermediate/filter_protein_groups/"
                                     "*_ProteinProphet_filtered.idXML")
                )
            inferencer = sorted(root.glob("*_IDConflictResolver.consensusXML"))
            if len(proteinprophet) + len(inferencer) != 1:
                raise ValueError("Expected one published DDA protein inference artifact")
            if proteinprophet:
                artifacts.append(self.artifact(proteinprophet[0], "idxml", "inference"))
            else:
                artifacts.append(self.artifact(inferencer[0], "consensusxml", "inference"))
            tmt_roots = sorted(tmt_directories[0].glob("OpenMS-*/output")) if tmt_directories else []
            if mode == "TMT" and not tmt_roots:
                warnings.append(
                    f"TMT quantitation directory {tmt_directories[0]} has no OpenMS output; "
                    "observation-level TMT quantities will be absent"
                )
            if tmt_roots:
                if len(tmt_roots) != 1:
                    raise ValueError(f"Expected one TMT output under {workflow}")
                mode_files = [(mode, sorted(tmt_roots[0].glob(f"*_{mode}.consensusXML")))
                              for mode in ("ScanBased", "FeatureBased")]
                found = [(mode, paths[0]) for mode, paths in mode_files if len(paths) == 1]
                if len(found) != 1 or any(len(paths) > 1 for _, paths in mode_files):
                    raise ValueError(f"Expected one scan or feature TMT consensusXML under {tmt_roots[0]}")
                quantification_mode, _ = found[0]
                resolved = sorted(tmt_roots[0].glob("*_IDConflictResolver.consensusXML"))
                if len(resolved) != 1:
                    raise ValueError(f"Expected one final TMT IDConflictResolver consensusXML under {tmt_roots[0]}")
                artifacts.append(self.artifact(
                    resolved[0], f"tmt_consensus_{'scan' if quantification_mode == 'ScanBased' else 'feature'}",
                    "quantitation"))
            artifacts.extend(self.qc_artifacts_dda(workflow, dataset))
        elif layer == "qc":
            artifacts.extend(self.qc_artifacts_dda(workflow, dataset))
        elif layer == "raw-search":
            warnings.append("DDA raw-search import requires run-specific idXML mapping")
        elif layer == "msbooster":
            warnings.append("DDA MSBooster import requires run-specific PIN mapping")
        if layer != "core" and not artifacts:
            raise ValueError(f"No supported files found for layer {layer} under {workflow}")
        return {**analysis, "layer": layer, "acquisition_mode": mode,
                "artifacts": artifacts, "warnings": warnings}

    def qc_artifacts_dda(self, workflow: Path, dataset_id: str) -> list[dict]:
        paths = sorted(workflow.glob(
            f"qc-ProteomicDDA-*/{dataset_id}/output/{dataset_id}_QC_metrics.tsv"
        ))
        if len(paths) > 1:
            raise ValueError(f"Multiple DDA QC tables under {workflow}: {paths}")
        return [self.artifact(path, "qc_metrics", "qc") for path in paths]

    def qc_artifacts(self, workflow: Path, dataset_id: str) -> list[dict]:
        """Select the single published cohort QC table, if available."""
        paths = sorted(
            workflow.glob(
                f"qc-ProteomicDIA-*/{dataset_id}/output/{dataset_id}_QC_metrics.tsv"
            )
        )
        if len(paths) > 1:
            raise ValueError(
                f"Multiple cohort QC tables found under {workflow}: {paths}"
            )
        return [self.artifact(path, "qc_metrics", "qc") for path in paths]

    @staticmethod
    def child_pipeline(workflow: Path, pattern: str) -> Path | None:
        matches = sorted(workflow.glob(pattern))
        if len(matches) > 1:
            raise ValueError(
                f"Multiple versions of {pattern} found; select one analysis"
            )
        return matches[0] if matches else None

    @staticmethod
    def artifact(path, format_name, stage, sample_name=None) -> dict:
        return dict(
            path=str(path), format=format_name, stage=stage, sample_name=sample_name
        )
