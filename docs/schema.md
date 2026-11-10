# ProteoVault schema 1.1.0

This schema stores one analysis per database. Metadata records its producer version, schema version, creation time, and import warnings. The reader validates the version before accepting it.

## Identity and provenance

Relationship IDs are stable analysis-local signed 63-bit BLAKE2b-derived integers over explicit entity fields. They are compact in joins and indexes. The raw source identifiers remain in their descriptive columns. They are not cross-analysis protein-group identifiers. PSM IDs identify a source file and its identification/hit ordinals, preserving multiple stages as observations. Source paths are also recorded relative to the supplied pipeline root. File size and modification time are recorded; they are not content checksums.

Unmodified sequence, modified sequence, and charged precursor are separate entities. Modified sequence identity includes its source notation. Protein accessions are preserved verbatim, including decoy/contaminant prefixes. Accession prefix flags are heuristic (`DECOY_`/`rev_` and text containing `contam`), not a replacement for source target/decoy annotations.

## Tables and views

- **`analysis`**: `analysis_id, dataset_id, pipeline, pipeline_version, root_path`.
- **`qc_metrics`**: source row, original metric dimensions, text value, unit, status, and source file from the merged cohort QC TSV.
- **`import_counts`**: `source_id, entity, row_count`.
- **`source_quant_methods`**: `source_id, level, method_id`; records available methods for each quantification source, including sources with zero quantity rows.
- **`metadata`**: `key, value`.
- **`peptide_proteins`**: `peptide_id, protein_id`.
- **`peptides`**: `peptide_id, sequence`.
- **`peptide_quant_features`**: one source ProteinQuantifier peptide row, including its exact modified `source_peptide`, reported protein, charge, and fraction. Its `peptide_id` links to the unmodified sequence in `peptides`.
- **`peptide_quant`**: `feature_id, channel_id, source_id, method_id, quantity`.
- **`peptide_quant_long`**: peptide quantities joined to sample names and source feature attributes.
- **`peptidoforms`**: `peptidoform_id, peptide_id, modified_sequence, notation`.
- **`pin_features`**: `source_id, row_number, spectrum_reference, scan_number, label, feature_name, value`.
- **`precursor_groups`**: `precursor_id, group_id`.
- **`precursor_quant`**: one row per source observation and channel, with `quant_id`, precursor/channel/source/method IDs, quantity, optional confidence, source element/map IDs, RT and peak boundaries in seconds, and an optional spectrum link.
- **`precursor_quant_long`**: readable precursor/sample/method context plus the observation ID, source element/map IDs, RT fields, and spectrum link.
- **`precursors`**: `precursor_id, peptidoform_id, charge`.
- **`protein_group_members`**: `group_id, protein_id`.
- **`protein_group_summary`**: `group_id, source_id, source_group_id, method, group_type, score, score_type, member_count`.
- **`protein_groups`**: `group_id, source_id, source_group_id, method_id, group_type, score, score_type`.
- **`protein_groups.reported_n_proteins` / `reported_n_peptides`**: source ProteinQuantifier row summaries when available.
- **`protein_quant`**: `group_id, channel_id, source_id, method_id, quantity`.
- **`protein_quant_long`**: `sample_name, group_id, source_group_id, grouping_method, quantification_method, quantity, source_id`.
- **`proteins`**: `protein_id, accession, sequence, description, is_decoy, is_contaminant`.
- **`protein_identifiers`**: `protein_id, source_id, source_accession, mapping_type, mapped_id`.
- **`psm_features`**: `psm_id, name, value`.
- **`psm_metadata`**: `psm_id, name, value`.
- **`psm_proteins`**: `psm_id, protein_id`.
- **`psm_scores`**: `psm_id, score_id, value, higher_better`.
- **`psms`**: `psm_id, spectrum_id, precursor_id, search_engine, candidate_rank, target_decoy, source_id, identification_index, hit_index, stage_id`.
- **`quant_channels`**: `channel_id, sample_id, run_id, label`.
- **`tmt_plexes`**: plex identity and reference channel.
- **`tmt_channels`**: reporter label and source description for a quant channel within a plex.
- **`tmt_runs`**: one plex assignment per existing acquisition run.
- **`runs`**: `run_id, sample_id, run_name, file_name`.
- **`sample_protein_counts`**: `sample_name, grouping_method, quantification_method, source_id, quantified_protein_groups`.
- **`samples`**: `sample_id, sample_name`.
- **`source_files`**: `source_id, path, relative_path, format, stage, size_bytes, mtime_ns, status, details`.
- **`spectra`**: `spectrum_id, run_id, native_id, rt_seconds, ms_level, precursor_mz, precursor_charge, source_id`.

## Row meanings and joins

- `psms` represents one source-stage candidate observation. Join `psm_scores`, `psm_features`, and `psm_metadata` by psm_id. Main scores carry higher_better. Numeric hit metadata remains in psm_features; original metadata values are retained as strings.
- `spectra` is unique by run_id/native_id. RT is in seconds for OpenMS idXML and mzML. DIA-NN report RT is not copied into this table because its scan/run association is not yet resolved.
- `peptide_proteins` records the union of identification-derived sequence mappings. `psm_proteins` preserves which proteins were attached to each individual PSM observation. DIA-NN quantitation does not create identification evidence or inference assignments.
- `protein_groups` records source-specific groups. Inference probability remains in score with its supplied score_type. OpenMS protein_group and indistinguishable_proteins definitions are kept as different group_type values. Quantitative groups use method DIA-NN and group_type quantification.
- `protein_group_members` permits multiple proteins per group. There is no implicit representative and no forced single-protein resolution.
- `precursor_groups` links DIA-NN precursors with a reported quantity to their quantitative groups. An identified peptide may have no quantitative group assignment. This link is not a claim about ProteinProphet peptide assignment.
- `quant_channels` is the measurement axis shared by protein, precursor, and peptide quantities. DIA uses LFQ channels. DDA LFQ uses one channel per sample; DDA TMT uses one per plex/reporter label, with TMT-only attributes in `tmt_channels`. TMT sample names are `<plex_id>|<reporter_label>`; channel descriptions are not used as sample identities. `tmt_runs` relates the plex to its mzML fractions, while TMT channels have no single `run_id`.
- `protein_quant` is unique by group/channel/source/method. The core imports `pg_matrix` as method `matrix` and available `PG.Quantity`, `PG.Normalised`, `PG.MaxLFQ`, and `PG.TopN` report columns under their original method names. Repeated precursor rows contribute only one value per group/run/method; conflicting nonmissing values reject the import. Matrix and report values remain distinct. One report source records both protein-group and precursor quantities.
- `protein_identifiers` stores one row per protein, source archive, mapping type, and mapped identifier. A UniProt idmapping archive may supply several values of the same type. The original `proteins.accession` remains unchanged. Only proteins with explicit `sp|accession|entry` or `tr|accession|entry` identifiers are matched; decoys, contaminants, and other identifier formats are skipped.
- `precursor_quant.quant_id` identifies one measured observation. DIA-NN matrix and report rows retain one value per precursor/channel/source/method; report `RT`, `RT.Start`, and `RT.Stop` are converted from minutes to seconds where available. Matrix rows have no RT. TMT scan-based rows retain one identified MS2 scan and reporter channel, with `spectrum_id` linked to its filtered PSM. Feature-based rows retain each consensus feature and source map; their centroid RT is stored, but `spectrum_id` is NULL because a feature can span several scans. `source_element_id` and `source_map_id` preserve repeated observations with the same precursor and channel. For each TMT source, `metadata` stores `tmt_quantity_scale:<source_id>` as `intensity`, `quasi_intensity`, or `log2_ratio`, determined from the consensusXML processing record. Intensities retain source units and should not be compared across methods without an explicit analysis choice.
- `pin_features` is deliberately source-row based. SpecId and ScanNr are retained without an unverified PSM link. Numeric PIN features can be queried immediately; a future exact bridge can link them when source identifiers support it.
- `import_counts` records expected imported PSM/quantity/feature counts per source. Validation compares these with database counts.

NULL means missing or unavailable. Zero is retained. No additional FDR threshold, decoy removal, contaminant removal, or quantification filtering is applied. Views include all imported rows. `sample_protein_counts` alone counts quantities greater than zero and groups by source and method. Always select a source/method explicitly before comparing matrices or inference results.

## Current limitations

No raw peak arrays, chromatograms, experimental sample annotations, cross-format peptidoform normalization, or arbitrary pipeline config parsing are imported automatically. DIA and DDA filtered sample idXML supply PSMs. DDA PSM runs are resolved from their `base_name` against the acquisition manifest; unmatched run names fail the build. Database identifiers must not be treated as equivalent across separately built analyses.

## Core and optional import layers

The schema remains the same across layers. A DIA core build populates DIA-NN protein and precursor quantities, validated PSMs, and protein inference groups. A DDA core build populates ProteinQuantifier protein and peptide quantities, filtered PSMs, inference groups, acquisition/run mappings, and available cohort QC. TMT builds also import identified reporter observations from the final `IDConflictResolver.consensusXML`, tagged as `TMT.ScanBased` or `TMT.FeatureBased`. Raw search PSMs and MSBooster PIN features remain optional for DIA. Core builds do not load full PSM feature/metadata annotations. Stage values in source_files and psms identify imported content; inspect derives the present layers from source provenance. Rebuild older databases to obtain the added tables.

FDR-controlled PSMs come from published `*_filter.idXML` files. ProteoVault preserves the pipeline's filtering, rather than applying its own threshold or substituting raw PSMs when filtered files are missing.

## Compact schema notes

`quant_methods`, `grouping_methods`, `processing_stages`, and `score_types` hold repeated labels once; large fact tables store their integer IDs. `score_aliases` records original idXML score names. Equal PEP aliases on the same PSM are stored as one score; differing values remain separate. Views ending in `_readable` expose text labels for queries, and the `*_long` views retain the R reader's interface. Databases built with the earlier layout need to be rebuilt; the reader supports only the compact layout.

Foreign-key and redundant natural-key indexes are omitted. `proteovault validate` checks relationships and natural-key uniqueness explicitly. Primary keys still guard fact and dimension deduplication. The integer key is a 63-bit digest; a theoretical hash collision is possible. The build validator detects duplicate natural keys and imported-record count mismatches, though it cannot prove the absence of every possible digest collision.

## Cohort QC metrics

If the published DIA or DDA cohort QC TSV exists during build, it is imported automatically. Otherwise, add it later with `proteovault add-layer DATABASE --layer qc`. The `qc_metrics` table preserves the source dimensions, including DDA `acquisition_mode`, `inference_route`, and `calculation_version`, plus `source_id` and source row number. The `value` column is text because some metrics are categorical. Empty values for unavailable metrics become SQL NULL; no QC values are recalculated or filtered.
