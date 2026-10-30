# Contributing

Install with `python -m pip install -e '.[dev-dependencies]'` and run `python -m pytest`.

Organize related operations into classes with coherent responsibilities. Keep each CLI's class, argument setup, and handler in its command module. Functions and methods are public by default; private names require a concrete reason. Use four-space indentation.

Format importers belong in separate source-format modules. Preserve source IDs, stages, shared evidence, missing values, and score meanings. Do not silently infer cross-format identifiers or change confidence filtering.

## Tests

Group tests into classes by domain. Place pure function and module tests in `test/unit/`. They should not read or create file fixtures. Place CLI and end-to-end database workflow tests in `test/integration/`. Keep representative input files in `test/fixtures/` and copy them to a temporary directory before a test changes them.
