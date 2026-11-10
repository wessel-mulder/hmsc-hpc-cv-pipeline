# Unusual temporal turnover simulations

## Independent-change null

The analysis compares recorded species presence/absence between the three Danish atlas periods on the same filtered set of grid cells and the full upstream species set. For each atlas pair, it conditions on the earlier occurrence matrix and the observed transition totals for every species.

For each species, the null randomly reallocates its observed losses among cells where it was recorded in the earlier atlas and reallocates its observed gains among cells where it was not recorded earlier. Losses and gains are sampled separately and without replacement. This preserves every species' earlier occupancy, number of losses, number of gains, and later occupancy. Species changes are randomized independently. The null does not constrain changes by geography, habitat, adjacency, or shared spatial patterns among species.

Cell TBI is `(losses + gains) / (species present in either atlas)`. Empty cells are assigned TBI zero. Each simulated site TBI is compared with its observed value, and counts strictly below, equal to, and above the observation are retained for later tail-probability calculations. Per-site null means and sample variances are stored. The simulation also stores each draw's grid-wide mean TBI and its summary statistics. It does not store the full site-by-simulation array.

The report tests both unusually high and unusually low TBI at every site. Inclusive lower and upper Monte Carlo tails include ties and use the add-one correction. The two-sided p-value is twice the smaller tail, capped at one; direction is determined by observed TBI minus null mean. The report shows a raw `p < 0.05` map for curiosity and a Holm–Bonferroni adjusted map across all site-by-atlas-pair tests. The unadjusted map does not account for multiplicity. The run stores these p-values, Holm adjusted values, and a direction/count summary in the RData file. `sitewise-turnover-unadjusted.png` and `sitewise-turnover-holm.png` show flagged sites in both directions: blue for less turnover than expected and orange for more; grey sites are not flagged. The color intensity shows the signed difference from the null mean. Each map includes a Bornholm inset, following the GWR mapping approach.

## Data and filtering

- Occurrences: `data-new/0-pre/raw/occurrences/Y_occurrences.csv`; retain all 201 species.
- Design: `data-new/0-pre/raw/design/studyDesign.csv`.
- Updated environmental land-use matrix: `data-new/0-pre/derived/environmental/X_environment.csv`. For each cell and atlas, sum `pct_class_0` (ocean) and `pct_class_77` (inland water); exclude cells above 75% in any atlas.
- Survey coverage: `data-new/0-pre/raw/effort/dækning atlas1+2.xlsx`; retain Good/Average coverage in Atlases 1 and 2. Atlas 3 coverage is treated as adequate because the workbook does not provide it.
- Use the intersection of qualifying cells across all three atlases. Missing or nonbinary occurrence values stop the run; they are not recoded as absences.

## Run settings and progress

The script defaults to 999,999 simulations per atlas pair, 10,000 draws per chunk, and six workers locally. On Slurm it uses `SLURM_CPUS_PER_TASK`; `--workers=N` overrides either setting. `--simulations=N`, `--chunk-size=N`, and `--seed=N` are also available for pilot runs or alternate reproducible runs. On non-Windows systems, chunks run with base R's `parallel::mclapply`; Windows runs serially.

Each completed chunk appends a row to `simulation-progress.log` in the active run folder. The log reports progress separately for Atlas pairs 1–2, 2–3, and 1–3, including completed draws, elapsed seconds, simulation rate, and estimated remaining seconds. With the default 10,000-draw chunk, it updates after each chunk completes for that pair.

## Results layout

Each run is saved in `simulation-results/<null-model>/<timestamp>_nsim-<count>_seed-<seed>/`. A run folder contains one `simulation-results.RData`, its progress log, a status marker, and that run's PNG figures. The site significance maps include a Bornholm inset. `null-distributions-all-pairs.png` shows simulated grid-wide mean TBI with a separate TBI x-axis for each atlas pair and a shared density y-axis. `observed-site-tbi-distributions.png` shows the observed site-level TBI distributions on a shared 0–1 x-axis and density scale; red lines mark the observed means. The RData stores run settings and input timestamps alongside the site-level and grid-wide simulation summaries. `simulation-results/latest-run.txt` points to the most recently completed run; the R Markdown reads this pointer and does not launch another simulation.

The current model folder is `simulation-results/independent-change/`. Future null models should use their own descriptive folder beneath `simulation-results/` so their outputs do not overwrite one another. Existing 9,999-draw results are retained in a legacy run folder and do not include the newer site-level tail counts.

Example:

```sh
Rscript scripts/1-atlas-analyses/unusual-temporal-turnover/unusual-temporal-turnover.R --workers=6 --simulations=999999 --chunk-size=10000 --seed=20261104
```
