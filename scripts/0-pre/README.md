# Current raw-to-model-input workflow

Run commands from the repository root. The current workflow uses only source
files under `data-new/0-pre/raw` and code under `scripts`; generated products
belong under `data-new/0-pre/derived` and `data-new/2-hmsc/derived`.

## One-command reproduction

```sh
# Fresh staging rebuild, without network access; compare to the saved outputs:
Rscript scripts/reproduce-preprocessing.R --verify --offline

# Rebuild the canonical derived outputs and current figures:
Rscript scripts/reproduce-preprocessing.R --offline
```

The 132 required CHELSAcruts monthly crops are already provided under raw.
Omit `--offline` to allow missing crops to be read and cropped directly from
the CHELSAcruts server. The offline rebuild is tested; downloading an entirely
new set of global files is not part of the reproducibility claim.

`--verify` starts in an empty staging folder: it does not reuse old derived
HILDA or climate tables. It rebuilds every current step, compares numerical
values and keys to the saved baseline, and records results, package versions
and file checksums in `data-new/reproducibility/`. Numeric comparisons use a
1e-7 tolerance for platform-dependent floating-point calculations. RData
compression bytes and PNG rendering are not required to be byte-identical.

## Individual steps and inputs

1. `preparing-traits.R`: raw `traits/birdTraits-merged-all.xlsx`,
   `sti/STI_Devictor.csv` and `taxonomy/IOC_Names_File_Plus-15.1_wide.csv` ->
   `Tr/Trait_table.csv`. Keep all source taxa at species level; collapse IOC
   subspecies to their parent and merge duplicate rows only when trait values
   agree. Record source rows, aliases and missingness. Corvus cornix receives
   the user-authorised Corvus corone STI proxy; no other related-taxon proxies.
2. `preparing-taxonomy.R`: that trait table and the raw IOC CSV ->
   `taxonomy/tree_fromPD.tre`, taxonomy crosswalk, distance matrix and PNG.
   All 226 species remain distinct; rank distances are not evolutionary times.
3. `preparing-environment.R`: raw DOF grid, HILDA annual state rasters and
   CHELSAcruts monthly crops -> `environmental/X_environment.csv`. All 2,165
   cells x 3 atlases are retained: 6,495 rows and 21 columns. No original HILDA
   classes are combined. `hh` is mean annual Simpson diversity over all 15
   classes, including ocean/freshwater. Percentages use all valid cell area.
   Only March-June temperature (degrees C) and precipitation (mm/season) are
   included. Temperature input /10; precipitation unscaled; -32768 is missing.
   HILDA periods: 1971-1974, 1993-1996, 2014-2017; climate Atlas 3: 2014-2016.
4. `../2-hmsc/S0_Data_Definitions.R`: new environmental and trait tables plus
   raw `occurrences/Y_occurrences.csv`, `design/studyDesign.csv` and the effort
   workbook -> the forest-split RData, filter audits and description.
5. `plot-environment.R`: current environmental CSV and raw grid ->
   `environmental/figures/current/`, with absolute/relative maps and a breeding
   temperature-precipitation point cloud. PNG only. Relative maps are z scores
   within atlas, calculated on the same complete cells across atlases.

Run individual scripts with `Rscript scripts/0-pre/<name>.R` (S0 is under
`scripts/2-hmsc`). Helpers are in `scripts/0-pre/environmental/`.
`preparing-environment.R --plan` lists its monthly inputs; `--refresh-landuse`
forces original-state extraction; otherwise its existing annual table may be
reused. The one-command verifier always forces the raw HILDA rebuild.

## S0 choices that are preserved explicitly

- Exclude site-atlas rows where ocean class 0 + freshwater class 77 exceeds
  75%; retain exactly 75%. Require complete values for the selected predictors.
- Good/Average coverage for Atlas 1 and corrected Atlas 2 `UNDS_KOD_1`.
  Atlas 3 retains the earlier assumed-Good rule; this is not measured coverage.
- Species require at least five occurrences in every atlas. No additional
  trait-completeness filter. No new requirement that all atlases share sites.
- Cropland = class 22 ONLY. Urban =11, pasture =33, grass/shrub =55.
  Deciduous forest =43+44; evergreen =42+41; mixed/other =40+45.
  Together with the two climate variables and heterogeneity (`hh`), these are ten model predictors.
  Cover predictors are percentage points (0-100); no scaling in S0.
- Heterogeneity (`hh`) is included unchanged from the environmental CSV. Classes
  23/24 are not separate cover predictors, but contribute to heterogeneity. Source CSVs are never rewritten by S0.
- S0 outputs X, Y, Tr, Design, with matching row/species orders. Coordinates
  remain projected EPSG:23032; the historical Design year labels are preserved.

Scientific unit checks:

```sh
Rscript scripts/0-pre/tests/check-environment.R
```

## Limits of this audit

The supplied occurrence and study-design CSVs are the starting inputs. Their
upstream compilation from original bird records was not recreated here.
Earlier expanded-variable figures in `environmental/figures/` are historical
explorations; `figures/current/` corresponds to the current final matrix.
S1 and subsequent model-fitting/prediction scripts have not been executed for
this revised dataset. S1 still needs its paths and predictor selection updated,
and a consistent scaling convention, before new models can be prepared.
