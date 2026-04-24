# Figure Folder

This folder separates figure code from rendered figure files.

## Folders

- `scripts/` contains scripts for main and supplementary manuscript figures.
- `exploratory/` contains figure experiments and model-comparison sketches.
- `outputs/main/` contains current main-figure candidates.
- `outputs/supplementary/` contains current supplementary-figure candidates.
- `outputs/exploratory/` contains additional rendered exploratory outputs.
- `outputs/assets/` contains helper image assets used by figure scripts.

No files were deleted during the publication cleanup. Outputs were grouped by
their apparent role based on file names (`fig`, `sfig`, `xfig`, and exploratory
prefixes).

Some scripts still save to paths like `misc-figures/<filename>.png`. If you
rerun them, they may recreate files at the top level unless their `ggsave()`
paths are updated.
