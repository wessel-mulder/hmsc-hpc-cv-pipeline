# Environmental variance partitioning through the three atlases

## Reference and adaptation

Antão et al. (2022), *Climate change reshuffles northern species within their niches*, Nature Climate Change 12:587-592, Figure 4 and Methods (PDF pages 5 and 8). DOI: https://doi.org/10.1038/s41558-022-01381-x.
The paper regresses species-level variance components, weighted by Tjur’s R², against decade within taxon and zone. Its public extraction script explicitly multiplies each species VP by that species Tjur’s R² (and rounds the result); this implementation retains full precision. Source: https://github.com/benweigel/RECChange/blob/main/03%20Model%20parameter%20extracts.R. That extraction script does not contain the Figure 4 trend fits; the regression specification is drawn from the paper.

## Project data and response

Uses the three 2026-07-27 Atlas1/2/3 MinOccs5 CoverageGoodAverage model outputs under `intermediate-outputs/hmsc-output/`. All 157 species are present in all three fitted species lists, with identical names and order after explicit matching. The eight environmental predictors are analysed individually; the site random effect is retained in diagnostic tables but excluded from environmental multiple-testing families.
For species s, predictor v and atlas t, W[s,v,t] = VP[s,v,t] * TjurR2[s,t]. The Tjur estimates are the stored full-model fitted discrimination values, not cross-validation scores. Each CSV VP value was verified against its saved VP.rds, and each Tjur value against MF.rds; component sums equal one. No models are refitted, and no predictors are rescaled.
W is a Tjur-weighted VP index, not a uniquely identified causal contribution or an exact decomposition of observed binary-response variance. Unweighted VP is a relative allocation among the model’s environmental and random components.

## Temporal model and uncertainty

Atlas order is coded 0, 1, 2 for the 1970s, 1990s and 2010s. Slopes are per atlas step (roughly twenty years), not per calendar decade. A slope is fitted within each species and averaged across the matched species. With these equally spaced three periods the slope is exactly (W2010s - W1970s)/2: the middle period does not determine the slope. The trajectories and all three paired contrasts are therefore also reported to expose reversals.
For this balanced panel, the average within-species slope equals the coefficient of pooled lm(W ~ atlas_step). The latter slope and its ordinary, unpaired p-value are exported for comparison with the paper-style regression; the ordinary p-value is not used for the plotted inference.
Primary intervals use 10,000 bootstrap resamples of whole species trajectories, with the same resampling indices across variables and periods (seed 20260909). Plotted q-values come from two-sided one-sample t-tests of species slopes against zero, Benjamini-Hochberg adjusted across eight environmental predictors. Relative-VP tests form a separate sensitivity family. Pairwise contrast tests are adjusted across 24 comparisons (eight predictors times three contrasts). The site random-effect tests are diagnostic and have no adjusted q-value.
Pairing handles repeated measurements of each species. The bootstrap and t-tests still treat species as independent/exchangeable replicates; they do not model phylogenetic dependence or dependence induced by a jointly fitted community model. They describe cross-species variation conditional on stored estimates and do not propagate HMSC posterior uncertainty, model-fit uncertainty or grid sampling uncertainty. Interpret significance as exploratory, not as a posterior probability of temporal change.

## Results

Atlas means below are W multiplied by 100; slopes and intervals are percentage points per atlas step.

| Predictor | 1970s | 1990s | 2010s | Slope per atlas step [95% bootstrap interval] | BH q |
|---|---:|---:|---:|---:|---:|
| Temperature | 1.277 | 1.612 | 1.129 | -0.074 [-0.198, +0.044] | 0.2391 |
| Precipitation | 1.096 | 0.719 | 1.257 | +0.080 [-0.054, +0.210] | 0.2391 |
| Land-use heterogeneity | 0.650 | 0.762 | 0.740 | +0.045 [-0.022, +0.108] | 0.2391 |
| Urban cover | 0.838 | 0.905 | 1.305 | +0.234 [+0.108, +0.378] | 0.002 |
| Cropland cover | 1.351 | 1.522 | 1.644 | +0.147 [+0.028, +0.265] | 0.02844 |
| Pasture cover | 0.252 | 0.399 | 0.623 | +0.186 [+0.131, +0.246] | 3.042e-08 |
| Forest cover | 1.713 | 2.286 | 2.651 | +0.469 [+0.278, +0.681] | 3.418e-05 |
| Grass / shrub cover | 0.260 | 0.362 | 0.406 | +0.073 [+0.035, +0.114] | 0.001052 |

Urban, cropland, pasture, forest and grass/shrub cover show positive weighted slopes with q < 0.05. Forest has the largest mean increase. Temperature, precipitation and land-use heterogeneity lack a clear weighted overall slope at this threshold. The same five cover variables increase in unweighted VP; heterogeneity also increases in unweighted VP. Thus the main cover-variable result is not generated solely by changing Tjur weights.

## Interpretation and checks

Increasing VP indicates greater allocation of fitted variation to that predictor under the model and environmental distribution of that atlas. It does not identify an increasing causal effect of land-use change. VP components share a denominator, environmental correlations can change, and the study-grid coverage differs across atlases. All these can affect comparisons. The spatial random component is included in the exports to inspect allocation changes.
The endpoint weighting decomposition is exact: change(VP * R2) = mean(endpoint R2) * change(VP) + mean(endpoint VP) * change(R2), calculated for each species and then averaged. This is an algebraic diagnostic, not a causal decomposition. The script checks that contributions sum to each weighted change, and that changes summed over all environmental and spatial components equal the change in mean Tjur R².

## Files and reproduction

Run `Rscript scripts/4-data-visualisation/atlas-vp-temporal-change.R` from the repository root.
Figures (PNG only) are in `manuscript-outputs/figures/exploratory/vp-temporal-change/`: trajectories, weighted/unweighted slope comparison, and paired period-change heatmap.
Analysis tables in this folder: `species-atlas-vp-values.csv`, `species-trends.csv`, `overall-trends.csv`, `atlas-summary.csv`, `paired-period-changes.csv`, and `weighted-change-decomposition.csv`. VP values are stored as fractions; changes, slopes and their intervals use the corresponding unmultiplied scale. Multiply these effect estimates and intervals by 100 to recover plotted units. Counts, indices, p-values and q-values retain their usual units. `input-provenance.csv` records source paths and MD5 checksums; `session-info.txt` records the R environment.
