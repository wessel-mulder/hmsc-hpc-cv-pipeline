# Intermediate Outputs

This folder stores compact intermediate products needed to reproduce current manuscript outputs.

- `hmsc-output/`: compact HMSC summaries copied from the current 2026-07-27 model results.
- `gwr-output/`: cached GWR fit objects and local coefficient summaries for the richness/turnover GWR figure.
- `reproducibility-index.csv`: reverse-output provenance index linking outputs, scripts, inputs, and local-only requirements.

Bulky fitted HMSC models, raw HPC chain outputs, and large prediction arrays remain local-only under `HmscOutputs/` for now. They are required by scripts that need full fitted posteriors, such as CWM maps, methods grid-cell reconstruction, and full-vs-CV accuracy extraction.
