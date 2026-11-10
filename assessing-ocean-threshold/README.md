# Assessing the ocean threshold

Run from the repository root:

```sh
Rscript assessing-ocean-threshold/assess-ocean-threshold.R
```

The script reads `data-new/0-pre/derived/environmental/X_environment.csv` and the Atlas 1–2 effort workbook. Land coverage is calculated for each site-atlas row as `100 - pct_class_0 - pct_class_77`, where the two HILDA classes represent ocean and inland water. The output CSV preserves each selected cell's component coverages, and a second CSV reports how many cells remain at minimum-land thresholds from 0% to 50%.

Atlas 1 and 2 include rows coded Good (1) using `UNDS_KODE_` and corrected `UNDS_KOD_1`, respectively. Atlas 3 has no measured effort estimate, so all available rows are retained under the existing assumed-good convention; this is explicitly labeled in the figures and tables. Percentages use available rows with finite land coverage from 0 to 100.

The current preprocessing rule excludes ocean + freshwater above 75%, equivalent to retaining land coverage **at least 25%**. The figure and summaries retain equality at the threshold. The distribution is descriptive; it does not assess whether the cutoff is ecologically optimal.

Outputs:

- `figures/land-coverage-distribution.png`: overlaid land-cover distributions, with the current threshold.
- `figures/threshold-sensitivity.png`: retained-cell percentage over alternative thresholds.
- `land-coverage-summary.csv`: per-atlas descriptive counts and land-cover summaries.
- `threshold-retention-summary.csv`: per-atlas counts/percentages for each threshold.
- `good-coverage-cell-land-cover.csv`: cell-level coverage values used in the plots.
