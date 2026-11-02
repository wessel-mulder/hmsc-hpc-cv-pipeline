# Chapter 1 HMSC Reproducibility Repository

This repository is being organized as the publication-facing reproducibility archive for Chapter 1. Its practical target is to reproduce the current manuscript figures and tables from documented derived data products plus local-only bulky model artifacts.

## Folder map

- `scripts/`: publication-facing scripts, ordered from preprocessing through HMSC fitting and figure/table generation.
- `data-shared/`: portable S0 RData (including taxonomy tree) tracked for running S1 on the HPC.
- `data-new/0-pre/raw/`: source-like inputs that are small enough to keep locally here. Check source permissions before public release.
- `data-new/0-pre/derived/`: lightweight analysis-ready products used by the HMSC prep and figure/table scripts.
- `intermediate-outputs/`: compact model summaries, GWR summaries/caches, and the reproducibility index.
- `manuscript-outputs/`: final manuscript figures and tables promoted from the working analysis folders.

The older `Data/`, `HmscOutputs/`, `misc-figures/`, and `.codex-work/` folders are migration and provenance sources. New publication-facing work should avoid writing there.

## Quick start

Run commands from the repository root.

```sh
Rscript scripts/4-post/fig1-vp-effect-sizes.R
Rscript scripts/4-post/fig2-trait-influence-gamma-atlas-support.R
Rscript scripts/4-post/methods-analysis-grid-cells.R
Rscript scripts/4-post/sfig-full-vs-cv-accuracy.R
Rscript scripts/4-post/sfig-raptor-species-full-tjur-prevalence-trends.R
Rscript scripts/4-post/sfig-component-influence-sankey-bump.R
Rscript scripts/4-post/fig10-gwr-richness-turnover-materials/fig10-gwr-richness-turnover.R
Rscript scripts/4-post/methods-audit/make_methods_audit_csvs.R
python3 scripts/4-post/methods-audit/make_supplementary_methods_tables_docx.py
```

Some scripts still require local-only fitted model artifacts under `HmscOutputs/`, especially scripts that read full fitted posteriors or prediction arrays. These are documented in `intermediate-outputs/README.md` and `intermediate-outputs/reproducibility-index.csv`.

## Current reproducibility status

The current raw-input-to-S0 workflow has been rebuilt and verified: species-level
traits, IOC taxonomy/tree, original-class environmental data, current diagnostic
figures and the 10-predictor forest-split Hmsc input RData. Results and commands
are documented in `scripts/0-pre/README.md` and `data-new/reproducibility/`.

```sh
Rscript scripts/reproduce-preprocessing.R --verify --offline
```

This offline command requires the local downloaded HILDA states and CHELSAcruts
crop cache. These bulky inputs and Hmsc chains/checkpoints/prediction arrays are
excluded from Git; see `scripts/0-pre/README.md`. The complete local `Data/` and `data-new/` folders are excluded from Git. Only the portable model inputs under `data-shared/` are shared. Original occurrence/design compilation is outside this
audit. S1 model construction has been checked for all three atlases using only the shared RData. Sampling and prediction have not been run; scaling still needs review.

## HPC model definitions

After pulling this repository on the HPC, run from the repository root:

```sh
Rscript scripts/2-hmsc/S1_Model_Definitions.R
```

S1 reads its RData, including the bundled tree, from `data-shared/`; no local `data-new/` files are
needed for this step. S0 refreshes these portable inputs when rerun locally.
S1 constructs unfitted models; sampling and job submission remain separate.
