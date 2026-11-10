# Methods Audit Tables

These scripts build supplementary methods tables that document grid-cell filtering, species filtering, and retained modelled species with trait/prevalence summaries.

Run from the repository root:

```sh
Rscript scripts/4-data-visualisation/methods-audit/make_methods_audit_csvs.R
python3 scripts/4-data-visualisation/methods-audit/make_supplementary_methods_tables_docx.py
```

Outputs are written to `manuscript-outputs/tables/methods-audit/`.
