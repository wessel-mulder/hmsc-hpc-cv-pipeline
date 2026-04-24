# Cleanup Notes

This repository was reorganized for publication readiness with conservative
moves only. No model outputs or figure files were deleted.

## Structural Changes

- Added `support_scripts/` for shared helper functions.
- Moved exploratory notebooks and rendered exploratory HTML files to
  `notebooks/`.
- Moved supervisor-meeting material to `docs/meeting-notes/`.
- Moved loose root-level figures to `outputs/figures/`.
- Reorganized `misc-figures/` into:
  - `scripts/`
  - `exploratory/`
  - `outputs/main/`
  - `outputs/supplementary/`
  - `outputs/exploratory/`
  - `outputs/assets/`

## Duplicated Helpers Centralized

- `support_scripts/data_helpers.R`
  - `read_effort_coverage()`
  - `clean_lulc_columns()`
  - `standardise_sti_species_names()`
- `support_scripts/project_paths.R`
  - `find_model_folders()`
  - `model_subdir()`
- `support_scripts/hmsc_helpers.R`
  - `fix_hpc_alpha_samples()`
  - `ecospat_boyce()`
- `support_scripts/plot_helpers.R`
  - `plot_model_fit_cv()`

## Remaining Good Cleanup Candidates

- Convert repeated `pattern2match` script settings into command-line arguments
  for all staged scripts.
- Decide whether `HmscOutputs/` should stay tracked in the publication repo or
  be archived externally.
- Check whether rendered HTML/PDF notebooks should be kept in git or regenerated
  from source.
- Update `misc-figures/scripts/` save paths if you want new figure renders to
  land directly in `misc-figures/outputs/`.
