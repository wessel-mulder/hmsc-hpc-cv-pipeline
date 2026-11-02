# Local inputs excluded from Git

The scripts, smaller raw inputs, derived CSV/RData files and current figures
are tracked. The following downloaded/generated bulky artifacts remain on the
working machine and are deliberately excluded from normal Git pushes:

- `0-pre/raw/hildap_vGLOB-2.0_geotiff_wgs84/`: HILDA+ v2 annual state/transition
  rasters. Obtain the annual state maps from https://doi.org/10.1594/PANGAEA.974335
  and place them under `states/` using their original filenames.
- `0-pre/raw/chelsacruts-v1-grid-crops/`: reusable cropped CHELSAcruts V1.0
  monthly source rasters. `scripts/0-pre/preparing-environment.R` populates
  missing crops directly from the source when run without `--offline`.
- Outside `data-new`, the `Models/`, `Tests/` and `Results/Preds/` directories
  beneath `intermediate-outputs/hmsc-output/`: model chains, checkpoints and
  prediction arrays. Compact summaries and figures remain tracked.

The offline reproduction check requires the local HILDA states and climate
crops. A new clone must obtain the HILDA states and populate the climate cache
before using `--offline`. The rebuild audit on the working machine verified
these files; the audit is not a claim that all large inputs are bundled in Git.
See `scripts/0-pre/README.md` for exact commands and temporal coverage.
