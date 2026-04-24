# Support Scripts

Shared helper functions live here so the staged scripts and exploratory
notebooks do not each need local copies.

## Files

- `data_helpers.R`: survey-effort loading, land-use column cleaning, and STI
  taxonomy harmonisation.
- `project_paths.R`: helpers for locating timestamped HMSC output folders.
- `hmsc_helpers.R`: HMSC postprocessing helpers, including HPC alpha conversion
  and Boyce index calculation.
- `plot_helpers.R`: common plotting helpers.

Load helpers with:

```r
source(file.path("support_scripts", "data_helpers.R"))
```

Scripts are written assuming the working directory is the repository root.
