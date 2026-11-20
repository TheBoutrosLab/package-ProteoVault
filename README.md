# ProteoVault

ProteoVault stores linked proteomics results in a portable DuckDB database for analysis in R, Python, and SQL. It imports completed metapipeline-ProteomicDIA analyses and DDA ProteinQuantifier LFQ or TMT abundance tables.

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

ProteoVault is licensed under the BSD 3-Clause License (`BSD-3-Clause`).
See [LICENSE](LICENSE) for the terms. Third-party components retain their respective licenses.

Copyright (C) 2026 Sanford Burnham Prebys Medical Discovery Institute ("Boutros Lab")
