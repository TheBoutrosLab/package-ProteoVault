# Interpretation and current coverage

- Missing quantities remain NULL; zero remains zero. Counts of quantified groups use abundance greater than zero.
- No new confidence filtering, normalization, or imputation is performed. Existing matrix filters remain those applied by DIA-NN; long reports retain supplied Q-values and PEPs.
- Protein groups from DIA-NN and inference idXML remain separate. A matching accession list does not establish identical inference provenance.
- The core reads PSMs from pipeline-filtered `*_filter.idXML` files and retains their main scores and recognized q-value/PEP annotations. It trusts the pipeline FDR filtering and does not impose another threshold. Full PSM feature/metadata tables are omitted from core imports. After adding `raw-search`, raw and validated PSMs remain separate observations linked to the same spectrum where identifiers match.
- OpenMS and DIA-NN modified sequences retain their source notation. Cross-format peptidoform equivalence is not inferred.
- PIN features preserve SpecId and ScanNr in `pin_features`. They are not guessed onto PSMs, because those identifiers may be rewritten by pipeline stages.
- Cohort inference idXML contributes proteins and groups. Its merged PSMs are omitted; sample-level idXML supplies PSM evidence.
- A published `manifests/quantitation.csv` maps acquisition filenames to canonical sample names; without it, quantification uses acquisition stems. Sample-level QC rows must match a quantified sample, and unmatched names fail import. A sample-level idXML is linked by an exact or unique underscore-bounded match; ambiguous matches fail. This initial DIA implementation models one sample per acquisition. TMT channel metadata and multifraction sample mapping need dedicated support.
- The optional `--fasta` build input adds full sequences for proteins found in the DIA analysis. mzML metadata remains an explicit importer and is not discovered automatically.
- DDA/TMT build discovery, consensusXML, OpenMS quantification CSV, standalone Fido/protXML, peptide-level summarization, and the R convenience package are future work.

See [the schema reference](schema.md) for table meanings and join rules.

## Schema version

Schema 1.0.0 uses compact integer relationship IDs and lookup tables for repeated labels. Existing databases built with the earlier layout must be rebuilt; schema 1.0.0 now denotes the compact layout. Direct SQL users can query `psms_readable`, `psm_scores_readable`, `protein_groups_readable`, `protein_quant_readable`, and `precursor_quant_readable` for text labels.
