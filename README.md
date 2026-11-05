# ProteoVault

ProteoVault stores linked proteomics results in a portable DuckDB database for analysis in R, Python, and SQL. It currently supports completed metapipeline-ProteomicDIA analyses.

Core builds retain available DIA-NN protein-group and precursor quantities from both the matrices and the detailed report.

Use `proteovault update-paths` to repair recorded source paths after moving a completed analysis; see the command-line guide for options and validation behavior.

## Documentation

Start with the [documentation overview](docs/index.md), then see:

- [Installation](docs/installation.md)
- [Container usage](docs/container.md)
- [Command-line guide](docs/cli.md)
- [Python interface](docs/python.md)
- [R and SQL access](docs/analysis.md)
- [Interpretation and current coverage](docs/interpretation.md)
- [Database schema](docs/schema.md)

## License

GNU GPL version 2 or later. Copyright 2026 Boutros Lab, Sanford Burnham Prebys Medical Discovery Institute.
