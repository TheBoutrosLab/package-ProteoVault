# ProteoVault documentation

ProteoVault links quantification, validated PSMs, protein groups, and source provenance in a portable DuckDB database. It currently imports completed metapipeline-ProteomicDIA outputs.

The default build imports DIA-NN protein-group and precursor quantification, FDR-filtered sample PSMs, and cohort protein-inference groups. The full DIA-NN report, raw search results, and MSBooster features are optional layers. A published merged cohort QC metrics table is imported automatically when present, or later with `add-layer`. Original source paths, pipeline stages, file sizes, and modification times are retained.

- [Installation and local development](installation.md)
- [Container usage](container.md)
- [Command-line guide](cli.md)
- [Python interface](python.md)
- [R and SQL access](analysis.md)
- [Interpretation and current coverage](interpretation.md)
- [Database schema](schema.md)
