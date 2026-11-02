# Chapter 1 Overview

This chapter uses repeated Danish breeding bird atlas data to ask how species
distributions and their environmental associations changed across three atlas
periods. The analysis is built around hierarchical modelling of species
communities, with the chapter acting as the modelling foundation for the broader
thesis work on bird redistribution, community change, and the mechanisms that
may structure those changes.

## Core Question

The chapter examines whether Danish breeding bird occurrence patterns are
associated with spatial structure, climate, and land use in consistent or
changing ways across atlas periods. It also asks whether species traits help
explain differences among species in their estimated responses.

In practical terms, the analysis compares:

- Atlas 1, representing the 1970s
- Atlas 2, representing the 1990s
- Atlas 3, representing the 2010s

The analysis is conducted on the Danish 5 km atlas grid after filtering to grid
cells with sufficient land area, environmental data, and survey coverage.

## Data Scope

The current analysis-ready HMSC input contains:

- 5,530 atlas-grid observations across the three atlas periods
- 1,776 retained grid cells in Atlas 1
- 1,850 retained grid cells in Atlas 2
- 1,904 retained grid cells in Atlas 3
- 157 retained breeding bird species
- 8 main modelled environmental predictors
- 3 trait variables used to structure species-level responses

Species filtering starts from 221 species, removes introduced or domestic taxa,
and excludes species with too few records for robust atlas-period modelling.

## Modelling Approach

The main analysis fits period-specific HMSC models to species occurrence
matrices. Each model uses a probit response for occurrence and includes:

- Climate predictors: breeding-season temperature and precipitation
- Land-use predictors: land-use heterogeneity plus urban, cropland, pasture, forest, and grass/shrub cover
- A spatial site-level random effect on the atlas grid
- Species traits: migration category, foraging guild, and species thermal index
- A phylogenetic structure derived from the local taxonomy workflow

The fitted models are used to estimate species-level environmental associations,
trait moderation of those associations, variance partitioning, predictive
accuracy, and community-level spatial outputs.

## Main Outputs

The current manuscript-facing outputs are organized around several linked figure sets:

- **Variance partitioning and effect sizes:** compares the relative contribution
  of spatial, climate, and land-use components across atlas periods.
- **Trait influence / Gamma support:** evaluates whether migration, foraging
  guild, and species thermal index help explain variation in species responses.
- **Spatial community outputs:** maps predicted richness and all-site beta
  turnover, then uses geographically weighted regression to summarize local
  associations with climate and land-use gradients.
- **Model performance diagnostics:** compares full-model and cross-validated
  accuracy, including AUC, Tjur R2, and RMSE.
- **Supplementary guild and prevalence checks:** explores whether performance
  and prevalence patterns differ for focal groups such as diurnal raptors.

The GWR workflow uses a fixed 100 km bisquare kernel as a pragmatic mesoscale
neighbourhood for summarizing local associations. This should be interpreted as
a scale for statistical smoothing and regional contrast, not as direct evidence
for a specific ecological mechanism.

## Repository Role

This repository is being shaped as a publication-facing reproducibility archive
for the chapter. Its current target is reproduction of the active manuscript
figures and tables from documented derived data products and compact model
summaries.

The main workflow is:

1. Prepare atlas, environmental, trait, taxonomy, and study-design inputs.
2. Define HMSC model objects for each atlas period.
3. Fit full and cross-validated HMSC models using the HPC workflow.
4. Extract model summaries, variance partitioning, accuracy metrics, and parameter estimates.
5. Generate manuscript PNG figures and supplementary tables.

Final manuscript outputs live in `manuscript-outputs/`, reusable compact
summaries live in `intermediate-outputs/`, and publication-facing scripts live
in `scripts/`.

## Current Reproducibility Caveat

The repository currently supports current-output reproducibility rather than a
complete raw-to-figure rebuild. Several upstream inputs are source-restricted,
absent, or local-only, including the original DOF atlas export, raw climate
rasters, HILDA land-use rasters, some trait inputs, and bulky fitted HMSC
objects.

Because of this, the clearest public-release framing is:

- code and small derived inputs in GitHub,
- compact summaries and final figures in the repository,
- bulky fitted models and prediction objects in an external archive,
- explicit data availability notes for source-restricted atlas and environmental inputs.
