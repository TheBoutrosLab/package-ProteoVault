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

proteovault add-layer /path/to/results.proteovault.duckdb --layer quant-details
proteovault add-layer /path/to/results.proteovault.duckdb --layer raw-search
proteovault add-layer /path/to/results.proteovault.duckdb --layer msbooster
proteovault add-layer /path/to/results.proteovault.duckdb --layer qc
proteovault add-fasta /path/to/results.proteovault.duckdb --fasta /path/to/search.fasta
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

The input directory must contain exactly one `metapipeline-ProteomicDIA-VERSION/DATASET/main_workflow` tree, or be the versioned directory itself. A published `quantitate-ProteomicDIA-*/DATASET/output/DIA-NN/*-report.pg_matrix.tsv` is required. Core also reads available protein-group quantities from `*-report.tsv`; detailed precursor quantities from that report remain optional. The core also discovers published `*_filter.idXML` sample files and supported cohort inference idXML. It never falls back to unfiltered search PSMs. Missing core artifacts are reported as warnings. Ambiguous pipeline versions are rejected.

`inspect` prints aligned tables by default; use `inspect --json` for structured output. Other commands return JSON on stdout. Timestamped INFO progress logs and warnings go to stderr, covering discovery, file imports with record counts, validation, and completion. Errors go to stderr with exit status 1; malformed CLI arguments use status 2. `--dry-run` creates no files. `--force` replaces an existing database only after the new database imports and validates successfully. Build on a local filesystem with atomic file publication support, then copy the closed database to shared storage if needed.

`add-fasta` adds protein sequences to an existing vault without rebuilding it.
It uses the same exact-accession matching as `build --fasta`, imports only
proteins already in the vault, and records the FASTA as a source. The command
copies the database to a temporary file, validates the updated copy, and then
replaces the original. Use `--dry-run` to check the number of proteins currently
missing sequences. An unchanged FASTA that was already imported is skipped;
a different reference FASTA requires a rebuild.

`add-layer` uses the pipeline output directory recorded in the database. Pass `--pipeline-output /new/path` after relocating it. Repeat `--layer` to add multiple layers together; `--dry-run` previews the files. Already imported, unchanged files are skipped. Changed previously imported files require a rebuild. Layers must belong to the same dataset and metapipeline version. The command imports into a temporary copy, validates it, and replaces the original only after success. `inspect` reports which layers are present. A QC table published after the core build can be added with `--layer qc`.

After moving the versioned result directory, run `update-paths` with its new location. It updates the analysis root, imported source paths, and QC provenance paths inside the old root. It also updates a FASTA stored inside the moved result tree; use `--fasta /new/path/search.fasta` when a separate FASTA moved. External input paths outside the result tree are left as recorded. Use `--dry-run` to preview source changes. The command checks source files and validates a temporary database copy before replacing the original.

| Layer | Contents |
| --- | --- |
| `core` (default build) | Quantification matrices, available `PG.Quantity`, `PG.Normalised`, `PG.MaxLFQ`, and `PG.TopN` report quantities, FDR-filtered PSMs with q-values/PEPs, inference groups |
| `quant-details` | Detailed DIA-NN report precursor quantities: original, normalized, and translated |
| `raw-search` | Raw search idXML, including its PSM metadata and numeric features |
| `msbooster` | Edited MSBooster PIN features |
| `qc` | Final merged `<dataset>_QC_metrics.tsv` from qc-ProteomicDIA |
