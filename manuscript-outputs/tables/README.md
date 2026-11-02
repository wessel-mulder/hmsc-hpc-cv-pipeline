# Tables

Current table outputs are stored under `methods-audit/` and `supplement/`.

Regenerate the CSVs with:

```sh
Rscript scripts/4-data-visualisation/methods-audit/make_methods_audit_csvs.R
```

Regenerate the Word table document with:

```sh
python3 scripts/4-data-visualisation/methods-audit/make_supplementary_methods_tables_docx.py
```

Regenerate the supplementary raptor accuracy permutation tests with:

```sh
Rscript scripts/4-data-visualisation/sfig3-diurnal-raptors-cv-accuracy-tests.R
```

Regenerate the supplementary guild temporal-change accuracy summaries and raptor-vs-guild permutation tests with:

```sh
Rscript scripts/4-data-visualisation/sfig4-guild-cv-accuracy-temporal-change.R
```

Regenerate the focused Tjur R2 temporal-change tests, prevalence-adjusted diagnostics, and guild-specific prevalence-slope tests with:

```sh
Rscript scripts/4-data-visualisation/sfig4-tjur-r2-guild-prevalence-tests.R
```
