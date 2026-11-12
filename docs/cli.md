# Command-line guide

```bash
proteovault build \
    --pipeline-output /path/to/output/metapipeline-ProteomicDIA-UniProtRP_ProteinProphet \
    --output /path/to/results.proteovault.duckdb \
    --dry-run

proteovault build \
    --pipeline-output /path/to/output/metapipeline-ProteomicDIA-UniProtRP_ProteinProphet \
    --fasta /path/to/the-search-target-decoy.fasta \
    --output /path/to/results.proteovault.duckdb

proteovault build \
    --pipeline-output /path/to/output/metapipeline-ProteomicDDA \
    --tmt-channel-csv /path/to/channel.csv \
    --output /path/to/dda.proteovault.duckdb

proteovault add-layer /path/to/results.proteovault.duckdb --layer raw-search
proteovault add-layer /path/to/results.proteovault.duckdb --layer msbooster
proteovault add-layer /path/to/results.proteovault.duckdb --layer qc
proteovault add-fasta /path/to/results.proteovault.duckdb --fasta /path/to/search.fasta
proteovault add-uniprot-idmapping /path/to/results.proteovault.duckdb --mapping /path/to/UP000005640_9606.idmapping.gz
proteovault update-paths /path/to/results.proteovault.duckdb --pipeline-output /new/path/to/metapipeline-ProteomicDIA-VERSION

proteovault validate /path/to/results.proteovault.duckdb
proteovault inspect /path/to/results.proteovault.duckdb
```

`--fasta` is optional. Supply the exact FASTA used for the search to populate
`proteins.sequence`. The build matches the first token of each FASTA header to
existing protein accessions and imports only matching sequences. It accepts plain
FASTA or gzip-compressed FASTA and records the file in `source_files`. A FASTA
with no matching accessions fails the build, because this usually means that
its identifiers or reference database do not match the search results. Without
`--fasta`, full protein sequences are available only when supplied by idXML.

The input directory must contain exactly one `metapipeline-ProteomicDIA-VERSION/DATASET/main_workflow` or `metapipeline-ProteomicDDA-VERSION/DATASET/main_workflow` tree, or be the versioned directory itself. DIA requires a published DIA-NN protein-group matrix. DDA requires the published ProteinQuantifier protein and peptide tables, an acquisition manifest (`manifests/qc-mzml.csv` or `manifests/validation-mzml.csv`), and a protein inference artifact. DDA TMT also requires `--tmt-channel-csv` with `plex_id,channel,description` columns; the metapipeline does not publish that input CSV into its result tree. Ambiguous pipeline versions are rejected.

For TMT, each reporter channel becomes one quantified sample named `<plex_id>|<reporter_label>` (for example `plex_01|126`); the source description remains in `tmt_channels`. The declared plex type in `run.resolved.yaml` supplies the physical reporter order. Reporters omitted from the channel CSV are still registered with an empty description so OpenMS abundance columns retain their correct mapping. `tmt_runs` links each mzML acquisition in the selected manifest to its plex. DDA core imports filtered PSMs by matching their `base_name` to an acquisition run. TMT also imports identified reporter observations from the final `IDConflictResolver.consensusXML` into `precursor_quant`; scan-based observations must match an accepted PSM, while feature-based observations retain their own element and map identities.
For multi-plex ProteinQuantifier exports, protein abundance columns must identify the plex and reporter channel explicitly, as in `abundance|plex_01|ch1`. Numbered `abundance_sampleN[unknown]` protein columns are rejected because their channel order cannot be verified from the peptide CSV. Use OpenMS ProteinQuantifier's `-file_and_channel_level_output true` option when generating protein quantities, and confirm that its file/channel granularity matches the intended analysis.

When `manifests/quantitation.csv` is present, `sample_name` and exactly one of `path` or `mzML` map DIA-NN acquisitions to biological sample names throughout quantification and QC. Column order and unrelated columns do not matter; `fraction_id` is accepted for fractionated acquisitions. Unknown or ambiguous names fail the import. Sample-level QC names must match quantified samples, while empty cohort-level QC names remain empty. Rebuild an existing vault to correct sample identities already stored in it.

`inspect` prints aligned tables by default; use `inspect --json` for structured output. Other commands return JSON on stdout. Timestamped INFO progress logs and warnings go to stderr, covering discovery, file imports with record counts, validation, and completion. Errors go to stderr with exit status 1; malformed CLI arguments use status 2. `--dry-run` creates no files. `--force` replaces an existing database only after the new database imports and validates successfully. Build on a local filesystem with atomic file publication support, then copy the closed database to shared storage if needed.

`add-fasta` adds protein sequences to an existing vault without rebuilding it.
It uses the same exact-accession matching as `build --fasta`, imports only
proteins already in the vault, and records the FASTA as a source. The command
copies the database to a temporary file, validates the updated copy, and then
replaces the original. Use `--dry-run` to check the number of proteins currently
missing sequences. An unchanged FASTA that was already imported is skipped;
a different reference FASTA requires a rebuild.

Use `build --uniprot-idmapping /path/to/file.idmapping.gz` or
`add-uniprot-idmapping` to import a local UniProt three-column mapping archive.
Only proteins with explicit UniProt FASTA identifiers are matched, and every
mapping type and value is retained in `protein_identifiers` with source
provenance. The original protein accession is unchanged. `add-uniprot-idmapping`
validates a temporary database copy before replacing the original; an unchanged
archive already imported is skipped. The command requires schema 1.1.0.

`add-layer` uses the pipeline output directory recorded in the database. Pass `--pipeline-output /new/path` after relocating it. Repeat `--layer` to add multiple layers together; `--dry-run` previews the files. Already imported, unchanged files are skipped. Changed previously imported files require a rebuild. Layers must belong to the same dataset and metapipeline version. The command imports into a temporary copy, validates it, and replaces the original only after success. `inspect` reports which layers are present. A QC table published after the core build can be added with `--layer qc`.

After moving the versioned result directory, run `update-paths` with its new location. It updates the analysis root, imported source paths, and QC provenance paths inside the old root. It also updates a FASTA stored inside the moved result tree; use `--fasta /new/path/search.fasta` when a separate FASTA moved. External input paths outside the result tree are left as recorded. Use `--dry-run` to preview source changes. The command checks source files and validates a temporary database copy before replacing the original.

| Layer | Contents |
| --- | --- |
| `core` (default build) | Quantification matrices, available `PG.Quantity`, `PG.Normalised`, `PG.MaxLFQ`, `PG.TopN`, `Precursor.Quantity`, `Precursor.Normalised`, and `Precursor.Translated` report quantities, FDR-filtered PSMs with q-values/PEPs, inference groups |
| `raw-search` | Raw search idXML, including its PSM metadata and numeric features |
| `msbooster` | Edited MSBooster PIN features |
| `qc` | Final merged `<dataset>_QC_metrics.tsv` from qc-ProteomicDIA |
| `uniprot-idmapping` | Protein identifier cross-references from a local UniProt archive |
