"""Discover the core or optional layers of one published DIA analysis."""

import re
from pathlib import Path


class PipelineOutputDiscovery:
    """Resolve one analysis and select artifacts for an explicitly requested layer."""

    LAYERS = ("core", "quant-details", "raw-search", "msbooster", "qc")

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve()

    def analysis(self) -> tuple[Path, dict]:
        if not self.root.is_dir():
            raise FileNotFoundError(
                f"Pipeline output directory does not exist: {self.root}"
            )
        if re.fullmatch(r"metapipeline-ProteomicDIA-[0-9].*", self.root.name):
            candidates = [self.root]
        else:
            candidates = sorted(
                p for p in self.root.glob("metapipeline-ProteomicDIA-*") if p.is_dir()
            )
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one metapipeline-ProteomicDIA version under {self.root}; found {len(candidates)}"
            )
        version_root = candidates[0]
        workflows = sorted(version_root.glob("*/main_workflow"))
        if len(workflows) != 1:
            raise ValueError(
                f"Expected one dataset/main_workflow under {version_root}; found {len(workflows)}"
            )
        workflow = workflows[0]
        return workflow, dict(
            dataset_id=workflow.parent.name,
            pipeline="metapipeline-ProteomicDIA",
            pipeline_version=version_root.name.removeprefix(
                "metapipeline-ProteomicDIA-"
            ),
            root=str(self.root),
        )

    def discover(self, layer: str = "core") -> dict:
        if layer not in self.LAYERS:
            raise ValueError(f"Unknown layer: {layer}")
        workflow, analysis = self.analysis()
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
                ("*-report.tsv", "diann_pg_report"),
            ]:
                matches = sorted(quant.glob(pattern))
                if len(matches) > 1:
                    raise ValueError(f"Ambiguous {format_name}: {matches}")
                artifacts.extend(
                    self.artifact(p, format_name, "quantitation") for p in matches
                )
                if format_name == "diann_pg_report" and not matches:
                    warnings.append("No DIA-NN report protein quantities found")
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
        elif layer == "quant-details":
            quant_roots = sorted(
                workflow.glob("quantitate-ProteomicDIA-*/*/output/DIA-NN")
            )
            if len(quant_roots) != 1:
                raise ValueError(
                    f"Expected one published DIA-NN output under {workflow}"
                )
            artifacts.extend(
                self.artifact(path, "diann_report", "quant-details")
                for path in sorted(quant_roots[0].glob("*-report.tsv"))
            )
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
