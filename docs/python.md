# Python interface

```python
from proteovault import VaultBuilder, VaultFastaAdder, VaultInspector, VaultLayerAdder, VaultValidator

builder = VaultBuilder(
    pipeline_output="/path/to/metapipeline-output",
    output="/path/to/results.proteovault.duckdb",
    fasta="/path/to/the-search-target-decoy.fasta",
)
plan = builder.plan()
result = builder.build()

layers = VaultLayerAdder(
    database="/path/to/results.proteovault.duckdb",
    layers=["raw-search", "msbooster"],
)
layers.add()

VaultFastaAdder(
    database="/path/to/results.proteovault.duckdb",
    fasta="/path/to/search.fasta",
).add()

report = VaultValidator(database="/path/to/results.proteovault.duckdb").validate()
summary = VaultInspector(database="/path/to/results.proteovault.duckdb").summary()
```

Each command class lives in its command module. Classes and ordinary methods are public by default.
