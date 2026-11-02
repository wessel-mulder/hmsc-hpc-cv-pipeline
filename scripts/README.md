# Chapter 1 scripts

Run from the repository root. Current folder structure:

- `0-pre/`: traits, IOC taxonomy/tree, environmental preparation and diagnostic plots.
- `1-atlas-analyses/`: atlas-level descriptive analyses, where present.
- `2-hmsc/`: S0 input assembly, S1 model definition and Hmsc/HPC workflows.
- `3-gwr/`: GWR workflow, where present.
- `4-post/`: post-processing and figure scripts (includes historical input paths).
- `5-vis/`: visualisation scripts.

The reproducible current workflow through S0 is described in
`0-pre/README.md`. Rebuild it with:

```sh
Rscript scripts/reproduce-preprocessing.R --offline
```

Validate a fresh staging rebuild against saved outputs with:

```sh
Rscript scripts/reproduce-preprocessing.R --verify --offline
```

Raw starting inputs: `data-new/0-pre/raw/`. Preprocessing products:
`data-new/0-pre/derived/`. S0 model inputs: `data-new/2-hmsc/derived/`.
Verification reports and checksums: `data-new/reproducibility/`.
This does not claim the historical downstream Hmsc/GWR/manuscript scripts are
ready for the revised inputs; see the reproduction report and S1 caveats.
