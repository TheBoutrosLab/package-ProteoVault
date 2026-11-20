# Changelog

All notable changes to ProteoVault are documented here.

## [Unreleased]

### Changed

- Migrate the package license from GPL to BSD 3-Clause and update license documentation and package metadata.

## [0.2.5] - 2026-11-19

## [0.2.4] - 2026-11-19

## [0.2.3] - 2026-11-13

### Fixed

- Import the 12-column QC v1.2.0 DDA metrics table directly for LFQ and TMT analyses.

## [0.2.2] - 2026-11-12

## [0.2.1] - 2026-11-10

## [0.2.0] - 2026-11-09

### Added

- Import DDA LFQ and TMT protein and peptide quantities, filtered PSMs, TMT reporter observations, and acquisition metadata into the core vault.

- Preserve DIA-NN report retention times and individual TMT scan or feature observations in `precursor_quant`.

- Import UniProt ID mapping archives into the provider-neutral `protein_identifiers` table during build or through `add-uniprot-idmapping` (schema 1.1.0).

- Accept published `sample_name,path` quantitation manifests and require sample-level QC names to resolve to quantified samples.

- Import distinct DIA-NN report protein-group quantities from `PG.Quantity`, `PG.Normalised`, `PG.MaxLFQ`, and `PG.TopN` in the core layer when present, collapsing repeated precursor rows with conflict checks.

- Add `update-paths` to repair recorded source and QC paths after moving a result directory.

## [0.1.0] - 2026-10-30

### Added

- `add-fasta` fills missing protein sequences in an existing vault from the search FASTA with atomic publication.

- Optional `build --fasta` imports matching protein sequences from the search FASTA and records the source file.

- Import the final merged cohort QC metrics TSV during build when present, or later with `add-layer --layer qc`.

- Core import of available DIA-NN report precursor quantities and confidence values.

- Compact schema 1.0.0 with integer relationship IDs and lookup tables for repeated score, method, and stage labels.

- Separate add-layer CLI and VaultLayerAdder for raw-search and MSBooster data, with dry runs, repeat-import skipping, and atomic publication.

- Inspection reports the layers present in a database.

- Timestamped CLI progress logging to stderr, preserving JSON-only stdout.

- Public VaultBuilder, VaultValidator, and VaultInspector classes with separate CLI modules.

- Atomic DuckDB builds with dry-run discovery and safe forced replacement.

- Versioned schema for sources, samples/runs, spectra, PSM observations and scores/features, peptide forms, protein groups, and quantities.

- DIA-NN matrix/report, OpenMS idXML, PIN, FASTA, and mzML metadata importers.

- Published metapipeline-ProteomicDIA output discovery and source provenance.

- Schema, relationship, quantity, and source-count validation; JSON inspection.

- R query examples and integration tests for import fidelity and failure handling.

### Changed

- `inspect` now prints readable tables by default; `--json` retains structured output.

- Core builds now import the DIA-NN matrices only; full report quantities are added explicitly.

- Remove foreign-key and redundant natural-key indexes; validate relationships and uniqueness during build.

- Collapse equal PEP score aliases while preserving original source names.

- Build now imports core quantification, FDR-filtered PSMs with confidence scores, and protein inference results; raw search and MSBooster feature data are opt-in.

- Cache repeated DIA-NN groups, precursors, channels, and mappings within bounded importer caches.

- Use typed Arrow schemas and larger batches to reduce Python conversion and database insert overhead.

### Fixed

- Allow pytest to import the local checkout without an editable package installation.

[0.1.0]: https://github.com/TheBoutrosLab/package-ProteoVault/releases/tag/v0.1.0
[0.2.0]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.1.0...v0.2.0
[0.2.1]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.0...v0.2.1
[0.2.2]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.1...v0.2.2
[0.2.3]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.2...v0.2.3
[0.2.4]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.3...v0.2.4
[0.2.5]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.4...v0.2.5
[unreleased]: https://github.com/TheBoutrosLab/package-ProteoVault/compare/v0.2.5...HEAD
