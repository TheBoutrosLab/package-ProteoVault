# Installation and local development

```bash
python -m pip install -e .
```

Python 3.10 or newer is required. DuckDB and PyArrow are installed as dependencies. For development:

```bash
python -m pip install -e '.[dev-dependencies]'
python -m pytest
```

The test suite groups related tests into domain-focused classes and separates file-free unit tests in `test/unit/` from CLI and database workflow tests in `test/integration/`. Prepared input files live in `test/fixtures/`; integration tests copy them to temporary directories before modification.
