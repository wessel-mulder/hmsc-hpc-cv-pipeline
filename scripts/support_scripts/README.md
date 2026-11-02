# Support Scripts

This folder contains small helper functions that are shared across multiple analysis or figure scripts.

## Contents

- `project_paths.R`: helpers for finding model folders and constructing model-output paths.

## Current Path Caveat

`project_paths.R` still defaults to the original `HmscOutputs/` layout. As the repository is migrated, this is the right place to centralize path changes so figure and fitting scripts can move toward `intermediate-outputs/`, `manuscript-outputs/`, and `data-new/` without each script hard-coding its own paths.

Keep this folder lightweight: shared path helpers, small formatting helpers, and reusable loaders belong here; figure-specific or model-specific code should stay beside the script that uses it.
