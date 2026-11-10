# S3 diagnostics validation, 2026-11-06

The complete combined workflow ran on all three `2026-11-02_06-55-21` atlas
fits, using all 250 saved draws in each of four chains and fitting thin = 100.
No fitted input was changed: input MD5 values matched before and after the run,
and were checked again against each report's provenance. Each model contains
155 species and 11 Beta coefficients (1,705 combinations).

| Model | All fitted cells | Shared cells | Flagged slopes / 1,550 | Flagged intercepts / 155 |
| --- | ---: | ---: | ---: | ---: |
| Atlas 1 | 1,747 | 1,697 | 60 | 144 |
| Atlas 2 | 1,821 | 1,697 | 568 | 147 |
| Atlas 3 | 1,875 | 1,697 | 684 | 155 |

No Beta combination had undefined diagnostics in these three fits. Passing is
reassuring diagnostic evidence, not proof of convergence. Fixed supplementary
parameters can have undefined diagnostics and are explicitly marked cannot assess.

## Checks completed

- Simulated independent, shifted-chain, autocorrelated, constant and invalid draws.
- Parameter order and species/predictor indexing checked against saved Beta draws.
- CSV median, SD, intervals, R-hat, bulk/tail ESS and MCSE (mean and endpoints)
  compared with direct `posterior` calls for three combinations per fitted model.
- Spatial low-rank correlations compared with direct `cor(X, Eta %*% Lambda)`
  for first/last draws and chains, on all and shared cells in each model.
- Combined spatial effects checked under compensating orthogonal rotations;
  projection and sum algebra checked independently.
- Environment-only, missing-trait, constant-predictor and undefined-diagnostic
  behavior exercised by a synthetic saved-object fixture through the complete
  legacy S3a entrypoint, including HTML/PDF rendering and validated MF cache reuse.
- Legacy `MF` object loading, cache fingerprint provenance, deterministic coefficient
  selections and shared-cell denominators checked against each actual model.
- HTML embeds PNGs and includes full-column table previews plus CSV links.
- PDFs knitted with XeLaTeX, rendered to PNG previews and visually inspected.
  Atlas PDFs have 93, 97 and 93 pages; the comparison has six pages. All atlas
  pages were rendered and reviewed as contact sheets, with selected pages checked
  at larger scale. Automated page checks found no blank pages or text outside
  page boundaries. PDF heatmaps were split to 32 species per page after inspection.
- Plot files are PNGs; PDFs are reports, explicitly requested by the user.
- R sources parse, and `git diff --check` passed. S1/S2 and fitted draws are unchanged.

## Reproduce

From the repository root:

```sh
CH1_MODEL_OUTPUT_DIR=data-new/2-hmsc/derived/models \
  Rscript scripts/2-hmsc/model-diagnostics/run-diagnostics.R 2026-11-02 250 100 4 4
Rscript scripts/2-hmsc/model-diagnostics/test-diagnostics.R --integration
CH1_MODEL_OUTPUT_DIR=data-new/2-hmsc/derived/models \
  Rscript scripts/2-hmsc/model-diagnostics/test-diagnostics.R --models
```

`Results/model-diagnostics/HPC_samples_0250_thin_100_chains_4/` beneath each
model contains the deliverables. The root's `model-diagnostics-comparison/`
contains the overview and combined directional-support family.

The initial prediction attempt failed because sandboxed CPU-core detection
returned NA. Running the same calculation outside the sandbox succeeded.
One PDF render failed during intermediate-file cleanup; staging a per-output
copy of the shared report template now keeps rendering intermediates isolated.
Pandoc emits a deprecated highlight-style option warning; the final report
renders and page checks pass. No statistical calculations were changed to work
around these rendering/environment issues.
