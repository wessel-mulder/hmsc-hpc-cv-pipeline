Species-environment Beta posterior contrasts across adjacent atlases
Model set: November 2, 2026; Atlas1/2/3 MinOccs5 CoverageGoodAverage.
Each model contains four chains x 250 retained draws; all retained draws were concatenated within atlas.
The preprocessed raw environmental matrix reproduces fitted XData after per-atlas standardization (max absolute error recorded in beta-scaling-audit.csv).
Each Beta draw was divided by the SD of its predictor among that atlas's exact fitted sites, yielding slopes in original predictor units.
Adjacent contrasts are Atlas2 minus Atlas1 and Atlas3 minus Atlas2.
For each contrast, 10000 posterior draw pairs were sampled independently from the two fitted atlas posteriors (seed 20261103).
PNG display filter: individual estimates with direction support above 85%. Full posterior tables and draw files are not filtered.
These are posterior contrasts between independently fitted models, not paired experimental tests.
Evidence is reported continuously using posterior medians, 95% credible intervals and sign probabilities; no binary significance cutoff or multiplicity adjustment is applied.
Only Betas with R-hat < 1.1 are included. Eligible species-predictor contrasts: 1513, 1429 of 1550 for Atlas1-2 and Atlas2-3.
beta-rhat-exclusion-audit.csv lists excluded species-predictor contrasts and the atlas model(s) that failed the threshold.
beta-change-community-wide-circles-support-over-85.png shows average slope changes across eligible species; only dots above the display support threshold are plotted.
beta-change-directional-species-examples-support-over-85.png names the three most positive and three most negative directional posterior contrasts per retained predictor and atlas pair.
Non-intercept Beta parameters with R-hat >= 1.1: 0, 37, 98 for Atlases 1, 2 and 3, respectively (out of 1,550 per atlas).
Treat contrasts involving poorly mixing parameters as exploratory; improve sampling/diagnostics before strong inferential claims.
Posterior uncertainty is conditional on the fitted model structures, predictor set, priors and training data.

Full paired posterior draws (local only): data-new/4-post/derived/compare-beta-changes/beta-adjacent-atlas-paired-posterior-draws.rds
