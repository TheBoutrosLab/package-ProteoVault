# R and SQL access

```r
connection <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = "/path/to/results.proteovault.duckdb",
    read_only = TRUE
)
protein_counts <- DBI::dbGetQuery(
    connection,
    "SELECT * FROM sample_protein_counts WHERE quantification_method = 'matrix'"
)
DBI::dbDisconnect(connection, shutdown = TRUE)
```

`protein_quant_long`, `precursor_quant_long`, `protein_group_summary`, and `sample_protein_counts` are ready-to-query views. Results retain source-specific protein groups and quantification methods; specify the relevant source/method when combining records.
