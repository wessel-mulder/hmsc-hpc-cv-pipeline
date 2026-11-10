# Exploring cells with supported excess temporal turnover

Generated: 2026-11-06 | Source run: 20261104T125100_nsim-999999_seed-20261104

## Main findings

All 30 flagged cells lost recorded species richness during their flagged interval. The clearest shared feature is an unusually large recorded loss, rather than a common extreme habitat profile. Richness and transition counts contribute directly to TBI, so these are descriptions of the selected outcome, not independent explanatory tests.

| Interval | Flagged cells | Median starting richness (others) | Median ending richness (others) | Median richness change (others) |
|---|---:|---:|---:|---:|
| 1-2 | 22 | 61.0 (53.0) | 33.0 (62.0) | -29.5 (8.0) |
| 2-3 | 6 | 71.0 (61.0) | 39.5 (63.5) | -27.5 (3.0) |
| 1-3 | 2 | 43.5 (53.0) | 25.5 (63.0) | -18.0 (10.0) |

All 22 Atlas 1→2 cells increase in richness again in Atlas 3. Nineteen finish Atlas 2 in the bottom tenth of the other cells’ richness distribution. Every flagged cell, in every interval, has losses in the top tenth and richness change in the bottom tenth of its interval’s reference distribution.

Among Atlas 1→2 flagged cells, the median fraction of lost species recorded again in Atlas 3 is 56.3%. Pooling species–cell losses gives 59.7% (484/811); these are different denominators.

### Survey coverage changes the interpretation

Nineteen of 22 Atlas 1→2 flagged cells (86%) have Average Atlas 2 coverage, compared with 484/1,676 other cells (29%). Seventeen of 22 (77%) shift from Good to Average, compared with 297/1,676 other cells (18%). The other three flagged cells retain Good coverage in both atlases. All six Atlas 2→3 flagged cells have Good Atlas 2 coverage; Atlas 3 effort is unavailable in this workbook.

The coverage association, richness rebound, and re-recording of lost species make uneven recording a concrete candidate explanation for the Atlas 2 trough. They do not establish that the changes are survey artefacts: real losses and subsequent recoveries could contribute. The current null does not model detection or effort.

### Habitat patterns

The 22 Atlas 1→2 cells span contrasting landscapes. Median starting total forest cover is 5.9% (other cells: 4.8%); median forest change is zero in both groups. Median annual-crop change is −5.2 percentage points (others: −4.3). These medians do not reveal a shared extreme forest/cropland trajectory; individual exceptions are retained in the heatmaps and cell profiles.

The six Atlas 2→3 cells are relatively forested: median starting forest cover is 25.7%, versus 5.7% elsewhere. Their median starting forest percentile is 81 nationally and 78 among the local references. Only one exceeds the national 90th percentile for total forest. Median forest change is −1.1 percentage points versus zero elsewhere. This is a small-group lead, not evidence that forest change explains their turnover.

The two Atlas 1→3 cells, BC45 and FG85, have contrasting trajectories and should be read individually. BC45 has 33→19→10 species; FG85 has 54→71→41.

Across all cells, the largest absolute descriptive Spearman association between excess turnover and the habitat variables considered here is about 0.15. This is weak evidence of a shared monotonic relationship; it does not rule out nonlinear, regional, or habitat-specific associations.

## Figures

1. [Richness trajectories](figures/01-richness-trajectories.png)
2. [Gains versus losses](figures/02-gains-versus-losses.png)
3. [Richness and transition percentiles](figures/03-richness-percentiles.png)
4. [Starting habitat percentiles](figures/04-starting-habitat-percentiles.png)
5. [Habitat change percentiles](figures/05-habitat-change-percentiles.png)
6. [Terrestrial-area cover sensitivity](figures/06-terrestrial-cover-percentiles.png)
7. [National versus local comparisons](figures/07-national-versus-local.png)
8. [Flagged cell locations, with Bornholm inset](figures/08-flagged-cell-locations.png)
9. [Continuous excess turnover](figures/09-continuous-excess-turnover.png)
10. [Richness change and survey coverage](figures/10-richness-and-survey-coverage.png)

## Definitions and reproducibility

Run from repository root: `Rscript scripts/1-atlas-analyses/unusual-temporal-turnover/explore-flagged-cells/explore-flagged-cells.R`.

The script loads the completed 999,999-draw simulation, checks all four current input modification times against its metadata, and uses its saved occurrence matrices and 1,698 common cells / 201 species. It checks reconstructed TBI and richness changes. It neither reruns nor modifies the simulation. Existing land/coverage inclusion rules are preserved: at least 25% terrestrial area in every atlas; Good/Average effort in Atlases 1 and 2; adequate Atlas 3 effort is an assumption.

Supported excess means positive observed TBI minus null mean and two-sided Holm-adjusted p < 0.05 over all 5,094 cell-by-interval tests. The independent-change null fixes each species’ earlier occupancy and gain/loss totals, but independently reallocates changes without geography, habitat, or effort constraints. No cells are flagged for excess in more than one interval.

TBI = (gains + losses) / species recorded in either atlas. Richness change = gains − losses. A cell’s baseline species composition influences its null expectation; excess is not raw turnover. “Other cells” includes all cells without supported excess, including cells with supported low turnover.

National reference percentiles compare each flagged cell with all other cells in its interval (1,676 / 1,692 / 1,696). Percentile = 100 × (number strictly below + half the number tied) / reference size. Bottom/top tenth use ≤10/≥90. Values are rounded to eight decimals for percentile and correlation ties to avoid numerical extraction noise. Zero-heavy cover variables can have broad tied ranks; an extreme percentile need not imply a large absolute change.

Whole-cell cover percentages include ocean and freshwater in the denominator. Total forest sums HILDA classes 40–45; annual crops means class 22 only. Terrestrial cover is 100 × period-mean cover / (100 − period-mean ocean − period-mean freshwater). Signed changes are later minus earlier, in percentage points. They measure net compositional change, not mapped gross habitat transitions. HILDA periods are 1971–1974, 1993–1996 and 2014–2017.

Existing hh is the mean annual Simpson diversity over all original classes including water. The additional terrestrial diversity is 1 − sum(p²) over terrestrial classes, calculated from period-mean cover renormalised to land. It is a descriptive sensitivity measure with a different temporal aggregation from hh. Neither is fragmentation.

Local references are exploratory comparisons, not a causal matching design or an exchangeability model. For each flagged cell and variable, consider other cells within 75 km. If fewer than 40 are available, use the nearest 100 other cells. Select 20 by smallest score: (distance/75 km)² plus squared differences in covariates divided by their whole-grid SDs. Richness/transition comparisons match on starting land area; habitat comparisons also match on starting richness; habitat changes additionally match on the same starting habitat variable. There are no strict covariate calipers, matching is variable-specific, and controls may be reused. Local land-area ranks are constrained by matching on land itself and should not be read as independent evidence. Exact controls, scores and distances are exported.

No new significance tests are performed. Spatial dependence, outcome-based selection, overlapping atlas intervals, imperfect matching and n=2/n=6 groups limit inference. Spearman correlations are descriptive; differences in selected cells do not independently validate the selection criterion.

Selection sensitivity compares Holm <0.05, raw p <0.05 and the top 5% of excess. Holm <0.01 selects zero because the smallest attainable initial Holm value with 999,999 draws and doubled inclusive Monte Carlo tails is 5,094 × 2/1,000,000 = 0.010188. This stricter cutoff is limited by Monte Carlo resolution, rather than demonstrating absence of stronger support.

## Exported data

- `data/all-cell-pair-profiles.csv`: all 5,094 rows, raw bird counts and every original cover class at start/end/change, plus aggregate/sensitivity variables and coverage.
- `data/flagged-cell-percentiles.csv`: flagged values, national/local percentiles and reference medians.
- `data/local-reference-cells.csv`: exact controls and distance/score for each variable and case.
- `data/group-descriptive-summary.csv`: medians and extreme-percentile counts by interval/variable.
- `data/flagged-cell-species-changes.csv`: each changed species in each flagged cell.
- `data/atlas2-losses-recorded-again.csv`: species lost in Atlas 2 and re-recorded in Atlas 3, for Atlas 1→2 flagged cells.
- `data/survey-coverage-summary.csv`, `data/continuous-associations.csv`, `data/selection-sensitivity.csv`.

## Individual cell profiles

Habitat entries below show the most extreme national reference percentiles among starting cover and cover changes (whole-cell metrics); percentile extremeness does not measure effect magnitude. Full profiles, including terrestrial denominators and local comparisons, are in the CSVs.

### BC43 | Atlas 1-2

- Richness: 65 → 28 → 59. Flagged interval: 2 gains, 39 losses, 26 persistent species.
- TBI 0.612; null mean 0.331; excess +0.281; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Freshwater: 19.37%; reference percentile 98.3; local percentile 95.0.
- Starting Needleleaf evergreen forest: 23.17%; reference percentile 96.1; local percentile 85.0.
- Starting Mixed forest: 4.64%; reference percentile 94.5; local percentile 95.0.
- Lost species: Anas platyrhynchos, Cygnus olor, Caprimulgus europaeus, Haematopus ostralegus, Thalasseus sandvicensis, Tringa totanus, Vanellus vanellus, Columba palumbus, Streptopelia decaocto, Falco tinnunculus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus schoenobaenus, Anthus trivialis, Coloeus monedula, Corvus cornix, Corvus corone, Curruca curruca, Cyanistes caeruleus, Delichon urbicum, Emberiza schoeniclus, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Lophophanes cristatus, Motacilla alba, Muscicapa striata, Oenanthe oenanthe, Passer montanus, Pica pica, Pyrrhula pyrrhula, Riparia riparia, Saxicola rubetra, Spinus spinus, Turdus viscivorus, Dendrocopos major, Picus viridis, Podiceps cristatus.
- Gained species: Cuculus canorus, Lanius collurio.

### BD26 | Atlas 1-2

- Richness: 65 → 33 → 67. Flagged interval: 4 gains, 36 losses, 29 persistent species.
- TBI 0.580; null mean 0.379; excess +0.200; Holm p 0.020256.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: -19.09 percentage points; reference percentile 1.6; local percentile 5.0.
- Change in Annual crops: 8.18 percentage points; reference percentile 96.9; local percentile 95.0.
- Starting Freshwater: 11.28%; reference percentile 94.6; local percentile 80.0.
- Lost species: Anas acuta, Aythya fuligula, Calidris pugnax, Chroicocephalus ridibundus, Gallinago gallinago, Larus canus, Streptopelia decaocto, Falco tinnunculus, Coturnix coturnix, Gallinula chloropus, Anthus trivialis, Carduelis carduelis, Coloeus monedula, Curruca curruca, Cyanistes caeruleus, Delichon urbicum, Hippolais icterina, Lophophanes cristatus, Lullula arborea, Luscinia luscinia, Motacilla alba, Motacilla flava, Muscicapa striata, Oenanthe oenanthe, Periparus ater, Phylloscopus collybita, Phylloscopus trochilus, Prunella modularis, Regulus regulus, Riparia riparia, Saxicola rubetra, Sylvia atricapilla, Sylvia borin, Troglodytes troglodytes, Turdus philomelos, Athene noctua.
- Gained species: Buteo buteo, Cygnus olor, Mergus serrator, Spatula querquedula.

### BF40 | Atlas 1-2

- Richness: 55 → 38 → 68. Flagged interval: 11 gains, 28 losses, 27 persistent species.
- TBI 0.591; null mean 0.365; excess +0.226; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: 13.09 percentage points; reference percentile 96.4; local percentile 90.0.
- Change in Annual crops: -14.55 percentage points; reference percentile 16.1; local percentile 30.0.
- Starting Pasture/rangeland: 5.62%; reference percentile 72.4; local percentile 55.0.
- Lost species: Gallinago gallinago, Gallinula chloropus, Acrocephalus palustris, Acrocephalus schoenobaenus, Anthus pratensis, Anthus trivialis, Carduelis carduelis, Coccothraustes coccothraustes, Curruca communis, Curruca curruca, Emberiza calandra, Emberiza schoeniclus, Hippolais icterina, Linaria cannabina, Locustella naevia, Lophophanes cristatus, Muscicapa striata, Periparus ater, Phylloscopus collybita, Phylloscopus sibilatrix, Prunella modularis, Sylvia atricapilla, Sylvia borin, Turdus viscivorus, Dendrocopos major, Asio otus, Athene noctua, Tyto alba.
- Gained species: Accipiter nisus, Buteo buteo, Falco tinnunculus, Fulica atra, Coloeus monedula, Corvus frugilegus, Lullula arborea, Phoenicurus ochruros, Phoenicurus phoenicurus, Spinus spinus, Strix aluco.

### BH59 | Atlas 1-2

- Richness: 55 → 30 → 59. Flagged interval: 10 gains, 35 losses, 20 persistent species.
- TBI 0.692; null mean 0.406; excess +0.286; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: -2.81 percentage points; reference percentile 13.6; local percentile 30.0.
- Starting Total forest: 0.00%; reference percentile 16.5; local percentile 15.0.
- Starting Unmanaged grass/shrub: 10.08%; reference percentile 79.4; local percentile 60.0.
- Lost species: Buteo buteo, Anas crecca, Gallinago gallinago, Haematopus ostralegus, Columba palumbus, Cuculus canorus, Gallinula chloropus, Acrocephalus schoenobaenus, Corvus corax, Cyanistes caeruleus, Delichon urbicum, Ficedula hypoleuca, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Lanius collurio, Lanius excubitor, Lophophanes cristatus, Loxia curvirostra, Motacilla alba, Oenanthe oenanthe, Parus major, Passer domesticus, Periparus ater, Phoenicurus phoenicurus, Pica pica, Prunella modularis, Regulus regulus, Saxicola rubetra, Sturnus vulgaris, Sylvia atricapilla, Sylvia borin, Turdus merula, Dendrocopos major, Asio otus.
- Gained species: Astur gentilis, Circus pygargus, Acanthis flammea, Acrocephalus palustris, Acrocephalus scirpaceus, Certhia brachydactyla, Emberiza schoeniclus, Muscicapa striata, Jynx torquilla, Tyto alba.

### CD89 | Atlas 1-2

- Richness: 62 → 34 → 66. Flagged interval: 7 gains, 35 losses, 27 persistent species.
- TBI 0.609; null mean 0.361; excess +0.248; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: 8.57 percentage points; reference percentile 92.3; local percentile 90.0.
- Change in Total forest: -1.57 percentage points; reference percentile 13.5; local percentile 15.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 19.9; local percentile 27.5.
- Lost species: Buteo buteo, Tadorna tadorna, Apus apus, Gallinago gallinago, Streptopelia decaocto, Alcedo atthis, Fulica atra, Anthus trivialis, Carduelis carduelis, Chloris chloris, Cinclus cinclus, Coloeus monedula, Corvus cornix, Corvus frugilegus, Curruca curruca, Cyanistes caeruleus, Garrulus glandarius, Hippolais icterina, Luscinia luscinia, Motacilla cinerea, Motacilla flava, Periparus ater, Phoenicurus phoenicurus, Pica pica, Poecile palustris, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Saxicola rubetra, Spinus spinus, Sylvia borin, Turdus viscivorus, Dendrocopos major, Asio otus, Strix aluco.
- Gained species: Acrocephalus schoenobaenus, Anthus pratensis, Curruca communis, Muscicapa striata, Phylloscopus collybita, Prunella modularis, Athene noctua.

### CE79 | Atlas 1-2

- Richness: 60 → 31 → 61. Flagged interval: 5 gains, 34 losses, 26 persistent species.
- TBI 0.600; null mean 0.296; excess +0.304; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Needleleaf evergreen forest: 1.00 percentage points; reference percentile 96.5; local percentile 100.0.
- Starting Urban: 11.00%; reference percentile 87.2; local percentile 95.0.
- Change in Pasture/rangeland: 5.50 percentage points; reference percentile 86.5; local percentile 75.0.
- Lost species: Buteo buteo, Anas platyrhynchos, Cygnus olor, Tadorna tadorna, Chroicocephalus ridibundus, Cuculus canorus, Gallinula chloropus, Acrocephalus scirpaceus, Aegithalos caudatus, Anthus trivialis, Carduelis carduelis, Coloeus monedula, Curruca communis, Curruca curruca, Delichon urbicum, Emberiza calandra, Emberiza schoeniclus, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Luscinia luscinia, Motacilla cinerea, Passer domesticus, Phylloscopus collybita, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Saxicola rubetra, Sitta europaea, Sylvia atricapilla, Sylvia borin, Dendrocopos major, Picus viridis, Tachybaptus ruficollis.
- Gained species: Perdix perdix, Cinclus cinclus, Lophophanes cristatus, Spinus spinus, Dryocopus martius.

### CE97 | Atlas 1-2

- Richness: 54 → 35 → 66. Flagged interval: 9 gains, 28 losses, 26 persistent species.
- TBI 0.587; null mean 0.355; excess +0.232; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Annual crops: 97.82%; reference percentile 87.9; local percentile 100.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 19.9; local percentile 17.5.
- Starting Pasture/rangeland: 0.00%; reference percentile 24.0; local percentile 25.0.
- Lost species: Cygnus olor, Tadorna tadorna, Apus apus, Cuculus canorus, Acrocephalus palustris, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Anthus trivialis, Curruca communis, Curruca curruca, Emberiza calandra, Emberiza schoeniclus, Garrulus glandarius, Hippolais icterina, Lanius collurio, Luscinia luscinia, Periparus ater, Poecile palustris, Prunella modularis, Pyrrhula pyrrhula, Saxicola rubetra, Sitta europaea, Sylvia atricapilla, Sylvia borin, Turdus merula, Tachybaptus ruficollis, Athene noctua, Strix aluco.
- Gained species: Buteo buteo, Perdix perdix, Certhia familiaris, Corvus corax, Ficedula hypoleuca, Phoenicurus ochruros, Phylloscopus sibilatrix, Dendrocopos major, Picus viridis.

### CG08 | Atlas 1-2

- Richness: 66 → 33 → 59. Flagged interval: 7 gains, 40 losses, 26 persistent species.
- TBI 0.644; null mean 0.348; excess +0.296; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Annual crops: -34.68 percentage points; reference percentile 1.0; local percentile 5.0.
- Change in Unmanaged grass/shrub: 22.30 percentage points; reference percentile 97.3; local percentile 95.0.
- Starting Annual crops: 100.00%; reference percentile 96.5; local percentile 97.5.
- Lost species: Accipiter nisus, Anas platyrhynchos, Tadorna tadorna, Gallinago gallinago, Numenius arquata, Vanellus vanellus, Streptopelia decaocto, Alcedo atthis, Fulica atra, Gallinula chloropus, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Anthus trivialis, Certhia familiaris, Corvus cornix, Corvus corone, Cyanistes caeruleus, Delichon urbicum, Emberiza schoeniclus, Garrulus glandarius, Linaria cannabina, Locustella naevia, Luscinia luscinia, Motacilla cinerea, Motacilla flava, Muscicapa striata, Passer domesticus, Periparus ater, Phylloscopus sibilatrix, Phylloscopus trochilus, Pica pica, Poecile palustris, Pyrrhula pyrrhula, Regulus regulus, Saxicola rubetra, Sitta europaea, Sturnus vulgaris, Ardea cinerea, Strix aluco, Tyto alba.
- Gained species: Haematopus ostralegus, Anthus pratensis, Certhia brachydactyla, Emberiza calandra, Phylloscopus collybita, Troglodytes troglodytes, Tachybaptus ruficollis.

### CG18 | Atlas 1-2

- Richness: 74 → 19 → 50. Flagged interval: 1 gains, 56 losses, 18 persistent species.
- TBI 0.760; null mean 0.329; excess +0.431; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: 18.70 percentage points; reference percentile 98.6; local percentile 100.0.
- Starting Annual crops: 100.00%; reference percentile 96.5; local percentile 95.0.
- Change in Annual crops: -20.01 percentage points; reference percentile 8.5; local percentile 30.0.
- Lost species: Anas platyrhynchos, Tadorna tadorna, Apus apus, Actitis hypoleucos, Gallinago gallinago, Tringa totanus, Vanellus vanellus, Ciconia ciconia, Columba palumbus, Streptopelia decaocto, Cuculus canorus, Falco tinnunculus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus scirpaceus, Carduelis carduelis, Certhia familiaris, Chloris chloris, Coccothraustes coccothraustes, Coloeus monedula, Corvus cornix, Corvus corone, Delichon urbicum, Emberiza schoeniclus, Erithacus rubecula, Galerida cristata, Garrulus glandarius, Hirundo rustica, Linaria cannabina, Locustella naevia, Lophophanes cristatus, Luscinia luscinia, Motacilla flava, Muscicapa striata, Oenanthe oenanthe, Passer domesticus, Passer montanus, Periparus ater, Phoenicurus phoenicurus, Phylloscopus collybita, Phylloscopus trochilus, Pica pica, Poecile palustris, Prunella modularis, Regulus regulus, Saxicola rubetra, Sitta europaea, Sturnus vulgaris, Troglodytes troglodytes, Turdus merula, Dendrocopos major, Podiceps cristatus, Asio otus, Athene noctua, Strix aluco.
- Gained species: Turdus viscivorus.

### CG27 | Atlas 1-2

- Richness: 75 → 38 → 61. Flagged interval: 5 gains, 42 losses, 33 persistent species.
- TBI 0.588; null mean 0.303; excess +0.284; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Broadleaf deciduous forest: -2.82 percentage points; reference percentile 0.7; local percentile 0.0.
- Change in Unmanaged grass/shrub: 23.33 percentage points; reference percentile 97.6; local percentile 85.0.
- Change in Annual crops: -25.44 percentage points; reference percentile 3.1; local percentile 30.0.
- Lost species: Pernis apivorus, Tadorna tadorna, Apus apus, Streptopelia decaocto, Alcedo atthis, Perdix perdix, Gallinula chloropus, Aegithalos caudatus, Carduelis carduelis, Certhia familiaris, Chloris chloris, Coccothraustes coccothraustes, Coloeus monedula, Corvus cornix, Curruca curruca, Delichon urbicum, Emberiza calandra, Emberiza schoeniclus, Ficedula hypoleuca, Garrulus glandarius, Hirundo rustica, Lanius collurio, Linaria cannabina, Locustella naevia, Lophophanes cristatus, Luscinia luscinia, Muscicapa striata, Oriolus oriolus, Passer domesticus, Periparus ater, Phylloscopus sibilatrix, Pica pica, Poecile palustris, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Sitta europaea, Spinus spinus, Turdus viscivorus, Dendrocopos major, Picus viridis, Asio otus.
- Gained species: Cygnus olor, Columba oenas, Certhia brachydactyla, Motacilla cinerea, Phoenicurus ochruros.

### CG28 | Atlas 1-2

- Richness: 65 → 26 → 56. Flagged interval: 3 gains, 42 losses, 23 persistent species.
- TBI 0.662; null mean 0.313; excess +0.349; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Annual crops: 100.00%; reference percentile 96.5; local percentile 95.0.
- Change in Pasture/rangeland: 9.59 percentage points; reference percentile 93.5; local percentile 90.0.
- Starting Total forest: 0.00%; reference percentile 16.5; local percentile 20.0.
- Lost species: Pernis apivorus, Anas crecca, Apus apus, Gallinago gallinago, Vanellus vanellus, Columba palumbus, Streptopelia decaocto, Alcedo atthis, Cuculus canorus, Falco tinnunculus, Fulica atra, Gallinula chloropus, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Carduelis carduelis, Certhia familiaris, Chloris chloris, Coccothraustes coccothraustes, Coloeus monedula, Corvus cornix, Curruca curruca, Cyanistes caeruleus, Delichon urbicum, Erithacus rubecula, Garrulus glandarius, Hirundo rustica, Linaria cannabina, Lophophanes cristatus, Luscinia luscinia, Motacilla alba, Muscicapa striata, Passer domesticus, Pica pica, Poecile palustris, Prunella modularis, Riparia riparia, Saxicola rubetra, Sitta europaea, Sturnus vulgaris, Turdus philomelos, Dendrocopos major, Strix aluco.
- Gained species: Certhia brachydactyla, Motacilla cinerea, Picus viridis.

### CH26 | Atlas 1-2

- Richness: 73 → 52 → 67. Flagged interval: 15 gains, 36 losses, 37 persistent species.
- TBI 0.580; null mean 0.315; excess +0.264; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Annual crops: -27.87 percentage points; reference percentile 2.3; local percentile 20.0.
- Change in Unmanaged grass/shrub: 23.04 percentage points; reference percentile 97.6; local percentile 90.0.
- Change in Pasture/rangeland: 4.82 percentage points; reference percentile 85.1; local percentile 80.0.
- Lost species: Circus aeruginosus, Circus pygargus, Anas crecca, Cygnus olor, Gallinago gallinago, Tringa totanus, Ciconia ciconia, Streptopelia decaocto, Fulica atra, Gallinula chloropus, Rallus aquaticus, Acrocephalus palustris, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Anthus pratensis, Anthus trivialis, Carduelis carduelis, Coloeus monedula, Curruca communis, Curruca curruca, Delichon urbicum, Emberiza schoeniclus, Locustella naevia, Luscinia luscinia, Motacilla flava, Phylloscopus sibilatrix, Phylloscopus trochilus, Poecile palustris, Regulus regulus, Riparia riparia, Saxicola rubetra, Sylvia borin, Dendrocopos major, Tachybaptus ruficollis, Asio flammeus, Asio otus.
- Gained species: Astur gentilis, Acanthis flammea, Certhia familiaris, Coccothraustes coccothraustes, Corvus corax, Ficedula hypoleuca, Lanius collurio, Loxia curvirostra, Phoenicurus ochruros, Poecile montanus, Pyrrhula pyrrhula, Regulus ignicapilla, Turdus pilaris, Ardea cinerea, Dryobates minor.

### CH66 | Atlas 1-2

- Richness: 72 → 50 → 58. Flagged interval: 11 gains, 33 losses, 39 persistent species.
- TBI 0.530; null mean 0.344; excess +0.186; Holm p 0.030354.
- Coverage: Good → Good → unavailable.
- Change in Total forest: -13.14 percentage points; reference percentile 1.8; local percentile 5.0.
- Change in Broadleaf deciduous forest: -1.24 percentage points; reference percentile 2.0; local percentile 10.0.
- Change in Annual crops: 9.37 percentage points; reference percentile 97.6; local percentile 85.0.
- Lost species: Aythya fuligula, Mergus serrator, Spatula querquedula, Apus apus, Chroicocephalus ridibundus, Gallinago gallinago, Haematopus ostralegus, Cuculus canorus, Perdix perdix, Rallus aquaticus, Acrocephalus scirpaceus, Coloeus monedula, Corvus cornix, Corvus corone, Curruca communis, Cyanistes caeruleus, Erithacus rubecula, Galerida cristata, Hippolais icterina, Locustella naevia, Luscinia luscinia, Motacilla cinerea, Motacilla flava, Oriolus oriolus, Parus major, Passer montanus, Phoenicurus phoenicurus, Poecile palustris, Pyrrhula pyrrhula, Riparia riparia, Sylvia borin, Troglodytes troglodytes, Picus viridis.
- Gained species: Astur gentilis, Pernis apivorus, Streptopelia decaocto, Falco tinnunculus, Corvus frugilegus, Curruca curruca, Emberiza calandra, Lanius collurio, Lophophanes cristatus, Phoenicurus ochruros, Strix aluco.

### CH95 | Atlas 1-2

- Richness: 77 → 37 → 61. Flagged interval: 0 gains, 40 losses, 37 persistent species.
- TBI 0.519; null mean 0.293; excess +0.226; Holm p 0.010188.
- Coverage: Good → Good → unavailable.
- Change in Annual crops: 3.08 percentage points; reference percentile 93.1; local percentile 75.0.
- Starting Pasture/rangeland: 19.77%; reference percentile 91.3; local percentile 90.0.
- Starting Annual crops: 3.79%; reference percentile 8.8; local percentile 5.0.
- Lost species: Accipiter nisus, Astur gentilis, Pernis apivorus, Mergus serrator, Tadorna tadorna, Charadrius hiaticula, Gallinago gallinago, Scolopax rusticola, Tringa totanus, Streptopelia decaocto, Falco tinnunculus, Perdix perdix, Acrocephalus palustris, Acrocephalus scirpaceus, Aegithalos caudatus, Anthus trivialis, Carduelis carduelis, Certhia brachydactyla, Certhia familiaris, Chloris chloris, Coloeus monedula, Corvus corax, Corvus cornix, Corvus corone, Corvus frugilegus, Curruca communis, Curruca curruca, Emberiza calandra, Ficedula hypoleuca, Hippolais icterina, Luscinia luscinia, Motacilla flava, Muscicapa striata, Oriolus oriolus, Phoenicurus phoenicurus, Phylloscopus sibilatrix, Poecile palustris, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus.
- Gained species: none.

### DF00 | Atlas 1-2

- Richness: 57 → 25 → 61. Flagged interval: 6 gains, 38 losses, 19 persistent species.
- TBI 0.698; null mean 0.378; excess +0.321; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Freshwater: 9.12%; reference percentile 92.7; local percentile 85.0.
- Change in Annual crops: 1.78 percentage points; reference percentile 90.7; local percentile 95.0.
- Change in Unmanaged grass/shrub: -1.78 percentage points; reference percentile 13.2; local percentile 15.0.
- Lost species: Cygnus olor, Mergus serrator, Tadorna tadorna, Charadrius hiaticula, Chroicocephalus ridibundus, Gallinago gallinago, Haematopus ostralegus, Tringa totanus, Vanellus vanellus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus schoenobaenus, Alauda arvensis, Anthus pratensis, Anthus trivialis, Carduelis carduelis, Corvus corax, Corvus cornix, Curruca curruca, Cyanistes caeruleus, Delichon urbicum, Emberiza calandra, Erithacus rubecula, Garrulus glandarius, Lophophanes cristatus, Motacilla flava, Passer montanus, Phoenicurus phoenicurus, Phylloscopus collybita, Prunella modularis, Riparia riparia, Sylvia atricapilla, Sylvia borin, Turdus philomelos, Dendrocopos major, Asio otus, Strix aluco.
- Gained species: Somateria mollissima, Falco tinnunculus, Acrocephalus scirpaceus, Coloeus monedula, Corvus frugilegus, Troglodytes troglodytes.

### DG24 | Atlas 1-2

- Richness: 52 → 37 → 80. Flagged interval: 11 gains, 26 losses, 26 persistent species.
- TBI 0.587; null mean 0.395; excess +0.192; Holm p 0.030354.
- Coverage: Good → Average → unavailable.
- Change in Broadleaf deciduous forest: 1.33 percentage points; reference percentile 95.8; local percentile 95.0.
- Starting Freshwater: 7.77%; reference percentile 89.6; local percentile 95.0.
- Change in Unmanaged grass/shrub: -1.33 percentage points; reference percentile 14.9; local percentile 5.0.
- Lost species: Anas crecca, Cygnus olor, Apus apus, Chroicocephalus ridibundus, Vanellus vanellus, Gallinula chloropus, Acrocephalus schoenobaenus, Anthus pratensis, Coloeus monedula, Corvus cornix, Emberiza schoeniclus, Erithacus rubecula, Ficedula hypoleuca, Hippolais icterina, Linaria cannabina, Phoenicurus phoenicurus, Phylloscopus sibilatrix, Pica pica, Poecile palustris, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Sitta europaea, Dendrocopos major, Tachybaptus ruficollis.
- Gained species: Accipiter nisus, Astur gentilis, Buteo buteo, Pernis apivorus, Aythya ferina, Aythya fuligula, Falco tinnunculus, Carduelis carduelis, Corvus corax, Periparus ater, Troglodytes troglodytes.

### DG59 | Atlas 1-2

- Richness: 60 → 45 → 56. Flagged interval: 12 gains, 27 losses, 33 persistent species.
- TBI 0.542; null mean 0.326; excess +0.215; Holm p 0.030354.
- Coverage: Good → Average → unavailable.
- Starting Urban: 10.05%; reference percentile 86.5; local percentile 75.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 19.9; local percentile 42.5.
- Starting Pasture/rangeland: 0.00%; reference percentile 24.0; local percentile 42.5.
- Lost species: Cygnus olor, Apus apus, Gallinago gallinago, Alcedo atthis, Cuculus canorus, Gallinula chloropus, Acrocephalus palustris, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Certhia familiaris, Curruca curruca, Emberiza calandra, Erithacus rubecula, Ficedula hypoleuca, Lophophanes cristatus, Luscinia luscinia, Phoenicurus phoenicurus, Phylloscopus sibilatrix, Pica pica, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Turdus philomelos, Dendrocopos major, Asio otus, Strix aluco.
- Gained species: Astur gentilis, Circus aeruginosus, Pernis apivorus, Tadorna tadorna, Scolopax rusticola, Tringa totanus, Carduelis carduelis, Certhia brachydactyla, Corvus corax, Periparus ater, Sitta europaea, Podiceps grisegena.

### FF08 | Atlas 1-2

- Richness: 55 → 36 → 44. Flagged interval: 13 gains, 32 losses, 23 persistent species.
- TBI 0.662; null mean 0.392; excess +0.269; Holm p 0.010188.
- Coverage: Good → Good → unavailable.
- Starting Mixed forest: 5.54%; reference percentile 95.7; local percentile 100.0.
- Starting Annual crops: 10.20%; reference percentile 11.8; local percentile 40.0.
- Starting Needleleaf evergreen forest: 1.66%; reference percentile 81.4; local percentile 80.0.
- Lost species: Cygnus olor, Somateria mollissima, Tadorna tadorna, Charadrius hiaticula, Haematopus ostralegus, Sterna hirundo, Tringa totanus, Vanellus vanellus, Falco tinnunculus, Perdix perdix, Gallinula chloropus, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Anthus pratensis, Carduelis carduelis, Chloris chloris, Coloeus monedula, Curruca communis, Curruca curruca, Delichon urbicum, Emberiza calandra, Emberiza schoeniclus, Linaria cannabina, Locustella naevia, Lullula arborea, Passer domesticus, Phoenicurus phoenicurus, Phylloscopus collybita, Phylloscopus trochilus, Prunella modularis, Riparia riparia, Sylvia atricapilla.
- Gained species: Cuculus canorus, Aegithalos caudatus, Carpodacus erythrinus, Corvus cornix, Corvus corone, Garrulus glandarius, Luscinia luscinia, Periparus ater, Pyrrhula pyrrhula, Regulus regulus, Sitta europaea, Spinus spinus, Sylvia borin.

### FF37 | Atlas 1-2

- Richness: 50 → 20 → 57. Flagged interval: 3 gains, 33 losses, 17 persistent species.
- TBI 0.679; null mean 0.407; excess +0.273; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Freshwater: 8.42%; reference percentile 91.2; local percentile 65.0.
- Change in Total forest: -3.28 percentage points; reference percentile 8.8; local percentile 10.0.
- Starting Broadleaf deciduous forest: 3.31%; reference percentile 88.5; local percentile 85.0.
- Lost species: Cygnus olor, Somateria mollissima, Tadorna tadorna, Gallinago gallinago, Haematopus ostralegus, Larus argentatus, Tringa totanus, Cuculus canorus, Gallinula chloropus, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Anthus trivialis, Carduelis carduelis, Curruca communis, Curruca curruca, Cyanistes caeruleus, Emberiza calandra, Emberiza schoeniclus, Erithacus rubecula, Hippolais icterina, Luscinia luscinia, Motacilla flava, Parus major, Passer montanus, Phoenicurus phoenicurus, Phylloscopus sibilatrix, Phylloscopus trochilus, Pica pica, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Sturnus vulgaris, Sylvia borin.
- Gained species: Streptopelia decaocto, Fulica atra, Delichon urbicum.

### FG71 | Atlas 1-2

- Richness: 58 → 15 → 47. Flagged interval: 0 gains, 43 losses, 15 persistent species.
- TBI 0.741; null mean 0.324; excess +0.417; Holm p 0.010188.
- Coverage: Average → Average → unavailable.
- Starting Broadleaf deciduous forest: 10.83%; reference percentile 96.2; local percentile 80.0.
- Change in Total forest: 5.33 percentage points; reference percentile 82.6; local percentile 70.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 19.9; local percentile 22.5.
- Lost species: Accipiter nisus, Buteo buteo, Cygnus olor, Scolopax rusticola, Columba palumbus, Streptopelia decaocto, Cuculus canorus, Falco tinnunculus, Perdix perdix, Gallinula chloropus, Acrocephalus palustris, Aegithalos caudatus, Carduelis carduelis, Chloris chloris, Delichon urbicum, Emberiza citrinella, Erithacus rubecula, Ficedula hypoleuca, Fringilla coelebs, Garrulus glandarius, Hippolais icterina, Lanius collurio, Luscinia luscinia, Motacilla alba, Passer montanus, Periparus ater, Phylloscopus collybita, Phylloscopus sibilatrix, Phylloscopus trochilus, Pica pica, Poecile palustris, Prunella modularis, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Sitta europaea, Sylvia atricapilla, Sylvia borin, Troglodytes troglodytes, Turdus philomelos, Dendrocopos major, Podiceps grisegena, Tachybaptus ruficollis.
- Gained species: none.

### FG81 | Atlas 1-2

- Richness: 62 → 13 → 52. Flagged interval: 1 gains, 50 losses, 12 persistent species.
- TBI 0.810; null mean 0.346; excess +0.464; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Starting Broadleaf deciduous forest: 13.55%; reference percentile 97.5; local percentile 90.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 19.9; local percentile 15.0.
- Starting Total forest: 17.58%; reference percentile 76.7; local percentile 75.0.
- Lost species: Buteo buteo, Spatula clypeata, Tadorna tadorna, Gallinago gallinago, Scolopax rusticola, Tringa totanus, Columba palumbus, Streptopelia decaocto, Cuculus canorus, Perdix perdix, Acrocephalus palustris, Acrocephalus schoenobaenus, Aegithalos caudatus, Anthus pratensis, Anthus trivialis, Carduelis carduelis, Certhia familiaris, Chloris chloris, Coloeus monedula, Curruca communis, Curruca curruca, Cyanistes caeruleus, Emberiza calandra, Emberiza citrinella, Emberiza schoeniclus, Erithacus rubecula, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Lanius collurio, Linaria cannabina, Luscinia luscinia, Motacilla flava, Oenanthe oenanthe, Parus major, Phoenicurus phoenicurus, Phylloscopus collybita, Phylloscopus sibilatrix, Phylloscopus trochilus, Pica pica, Poecile palustris, Prunella modularis, Regulus regulus, Riparia riparia, Sitta europaea, Sylvia atricapilla, Sylvia borin, Troglodytes troglodytes, Turdus philomelos, Dendrocopos major.
- Gained species: Anser anser.

### FG94 | Atlas 1-2

- Richness: 43 → 16 → 69. Flagged interval: 11 gains, 38 losses, 5 persistent species.
- TBI 0.907; null mean 0.442; excess +0.465; Holm p 0.010188.
- Coverage: Average → Average → unavailable.
- Change in Pasture/rangeland: -5.67 percentage points; reference percentile 8.8; local percentile 0.0.
- Change in Total forest: 7.29 percentage points; reference percentile 89.7; local percentile 70.0.
- Change in Unmanaged grass/shrub: -0.24 percentage points; reference percentile 18.1; local percentile 10.0.
- Lost species: Vanellus vanellus, Columba palumbus, Cuculus canorus, Alauda arvensis, Anthus trivialis, Carduelis carduelis, Chloris chloris, Corvus cornix, Curruca communis, Curruca curruca, Cyanistes caeruleus, Emberiza calandra, Emberiza citrinella, Erithacus rubecula, Fringilla coelebs, Hippolais icterina, Linaria cannabina, Luscinia luscinia, Motacilla alba, Muscicapa striata, Oriolus oriolus, Parus major, Periparus ater, Phylloscopus collybita, Phylloscopus sibilatrix, Phylloscopus trochilus, Poecile palustris, Prunella modularis, Regulus regulus, Riparia riparia, Saxicola rubetra, Sitta europaea, Sylvia atricapilla, Sylvia borin, Troglodytes troglodytes, Turdus philomelos, Dendrocopos major, Strix aluco.
- Gained species: Anas platyrhynchos, Mergus serrator, Tadorna tadorna, Streptopelia decaocto, Falco tinnunculus, Coloeus monedula, Corvus corax, Delichon urbicum, Locustella naevia, Passer montanus, Pica pica.

### BF05 | Atlas 2-3

- Richness: 59 → 91 → 66. Flagged interval: 12 gains, 37 losses, 54 persistent species.
- TBI 0.476; null mean 0.291; excess +0.184; Holm p 0.030354.
- Coverage: Good → Good → unavailable.
- Change in Needleleaf evergreen forest: 5.77 percentage points; reference percentile 96.6; local percentile 55.0.
- Starting Needleleaf evergreen forest: 22.38%; reference percentile 96.3; local percentile 85.0.
- Starting Unmanaged grass/shrub: 28.87%; reference percentile 94.3; local percentile 85.0.
- Lost species: Accipiter nisus, Anas platyrhynchos, Anser anser, Tadorna tadorna, Apus apus, Haematopus ostralegus, Tringa totanus, Streptopelia decaocto, Falco tinnunculus, Perdix perdix, Gallinula chloropus, Aegithalos caudatus, Anthus pratensis, Carduelis carduelis, Coloeus monedula, Corvus cornix, Corvus corone, Curruca curruca, Delichon urbicum, Emberiza schoeniclus, Ficedula hypoleuca, Garrulus glandarius, Hirundo rustica, Locustella naevia, Loxia curvirostra, Motacilla flava, Muscicapa striata, Oenanthe oenanthe, Passer domesticus, Poecile palustris, Pyrrhula pyrrhula, Riparia riparia, Sylvia borin, Turdus pilaris, Podiceps cristatus, Tachybaptus ruficollis, Asio otus.
- Gained species: Haliaeetus albicilla, Aythya fuligula, Numenius arquata, Alcedo atthis, Grus grus, Luscinia svecica, Poecile montanus, Remiz pendulinus, Saxicola rubicola, Botaurus stellaris, Jynx torquilla, Aegolius funereus.

### FG84 | Atlas 2-3

- Richness: 35 → 57 → 36. Flagged interval: 7 gains, 28 losses, 29 persistent species.
- TBI 0.547; null mean 0.336; excess +0.211; Holm p 0.030354.
- Coverage: Average → Good → unavailable.
- Change in Broadleaf deciduous forest: 1.55 percentage points; reference percentile 96.4; local percentile 80.0.
- Starting Broadleaf deciduous forest: 8.03%; reference percentile 93.5; local percentile 80.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 10.1; local percentile 15.0.
- Lost species: Astur gentilis, Buteo buteo, Anas platyrhynchos, Scolopax rusticola, Streptopelia decaocto, Falco tinnunculus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus scirpaceus, Aegithalos caudatus, Anthus trivialis, Certhia familiaris, Corvus cornix, Erithacus rubecula, Ficedula hypoleuca, Garrulus glandarius, Hirundo rustica, Motacilla alba, Periparus ater, Phylloscopus sibilatrix, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Sitta europaea, Sturnus vulgaris, Sylvia borin, Dendrocopos major.
- Gained species: Circus aeruginosus, Anthus pratensis, Coccothraustes coccothraustes, Corvus frugilegus, Delichon urbicum, Podiceps grisegena, Asio otus.

### FG97 | Atlas 2-3

- Richness: 51 → 77 → 42. Flagged interval: 3 gains, 38 losses, 39 persistent species.
- TBI 0.512; null mean 0.293; excess +0.220; Holm p 0.010188.
- Coverage: Average → Good → unavailable.
- Change in Total forest: -2.13 percentage points; reference percentile 11.1; local percentile 25.0.
- Starting Broadleaf deciduous forest: 3.48%; reference percentile 88.2; local percentile 50.0.
- Starting Pasture/rangeland: 0.00%; reference percentile 19.7; local percentile 27.5.
- Lost species: Accipiter nisus, Circus aeruginosus, Pernis apivorus, Anas platyrhynchos, Anser anser, Cygnus olor, Tadorna tadorna, Gallinago gallinago, Haematopus ostralegus, Scolopax rusticola, Tringa totanus, Columba oenas, Streptopelia decaocto, Falco tinnunculus, Rallus aquaticus, Aegithalos caudatus, Anthus pratensis, Anthus trivialis, Certhia familiaris, Coloeus monedula, Corvus frugilegus, Curruca curruca, Emberiza schoeniclus, Ficedula hypoleuca, Garrulus glandarius, Hippolais icterina, Loxia curvirostra, Muscicapa striata, Periparus ater, Phoenicurus ochruros, Pica pica, Pyrrhula pyrrhula, Regulus regulus, Saxicola rubetra, Sitta europaea, Sturnus vulgaris, Turdus viscivorus, Tachybaptus ruficollis.
- Gained species: Acrocephalus palustris, Emberiza calandra, Turdus pilaris.

### FH65 | Atlas 2-3

- Richness: 15 → 68 → 49. Flagged interval: 11 gains, 30 losses, 38 persistent species.
- TBI 0.519; null mean 0.314; excess +0.205; Holm p 0.030354.
- Coverage: Good → Good → unavailable.
- Starting Freshwater: 17.46%; reference percentile 97.9; local percentile 90.0.
- Change in Total forest: -5.72 percentage points; reference percentile 2.9; local percentile 5.0.
- Change in Urban: -0.95 percentage points; reference percentile 3.8; local percentile 0.0.
- Lost species: Accipiter nisus, Buteo buteo, Anas platyrhynchos, Cygnus olor, Mergus serrator, Spatula clypeata, Recurvirostra avosetta, Scolopax rusticola, Tringa totanus, Vanellus vanellus, Perdix perdix, Acrocephalus palustris, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Aegithalos caudatus, Anthus trivialis, Carduelis carduelis, Certhia familiaris, Corvus cornix, Corvus frugilegus, Curruca curruca, Emberiza schoeniclus, Ficedula hypoleuca, Garrulus glandarius, Hippolais icterina, Linaria cannabina, Luscinia luscinia, Phoenicurus ochruros, Pica pica, Strix aluco.
- Gained species: Milvus milvus, Pernis apivorus, Mergus merganser, Apus apus, Certhia brachydactyla, Corvus corax, Periparus ater, Phoenicurus phoenicurus, Regulus ignicapilla, Ardea cinerea, Tyto alba.

### GG05 | Atlas 2-3

- Richness: 61 → 64 → 34. Flagged interval: 2 gains, 32 losses, 32 persistent species.
- TBI 0.515; null mean 0.315; excess +0.200; Holm p 0.040424.
- Coverage: Good → Good → unavailable.
- Change in Total forest: -11.33 percentage points; reference percentile 0.7; local percentile 5.0.
- Starting Broadleaf deciduous forest: 16.70%; reference percentile 98.1; local percentile 100.0.
- Starting Total forest: 61.12%; reference percentile 97.0; local percentile 100.0.
- Lost species: Accipiter nisus, Circus aeruginosus, Anas platyrhynchos, Mergus merganser, Tadorna tadorna, Haematopus ostralegus, Scolopax rusticola, Columba oenas, Streptopelia decaocto, Falco tinnunculus, Perdix perdix, Gallinula chloropus, Acrocephalus palustris, Anthus trivialis, Certhia familiaris, Corvus corax, Corvus cornix, Corvus frugilegus, Ficedula hypoleuca, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Motacilla alba, Muscicapa striata, Passer domesticus, Phylloscopus sibilatrix, Phylloscopus trochilus, Pyrrhula pyrrhula, Regulus regulus, Riparia riparia, Ardea cinerea, Strix aluco.
- Gained species: Coccothraustes coccothraustes, Sylvia borin.

### GG15 | Atlas 2-3

- Richness: 55 → 74 → 37. Flagged interval: 4 gains, 41 losses, 33 persistent species.
- TBI 0.577; null mean 0.338; excess +0.239; Holm p 0.010188.
- Coverage: Good → Good → unavailable.
- Starting Broadleaf deciduous forest: 8.79%; reference percentile 94.4; local percentile 80.0.
- Starting Unmanaged grass/shrub: 0.00%; reference percentile 10.1; local percentile 10.0.
- Starting Freshwater: 5.67%; reference percentile 84.9; local percentile 60.0.
- Lost species: Accipiter nisus, Astur gentilis, Milvus milvus, Anas platyrhynchos, Cygnus olor, Mergus serrator, Somateria mollissima, Tadorna tadorna, Apus apus, Charadrius hiaticula, Gallinago gallinago, Haematopus ostralegus, Scolopax rusticola, Streptopelia decaocto, Cuculus canorus, Falco tinnunculus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus palustris, Aegithalos caudatus, Anthus trivialis, Certhia familiaris, Coccothraustes coccothraustes, Coloeus monedula, Curruca curruca, Delichon urbicum, Ficedula parva, Garrulus glandarius, Hippolais icterina, Hirundo rustica, Nucifraga caryocatactes, Oriolus oriolus, Periparus ater, Phoenicurus phoenicurus, Pica pica, Regulus regulus, Riparia riparia, Sylvia atricapilla, Sylvia borin, Turdus pilaris.
- Gained species: Haliaeetus albicilla, Anthus pratensis, Regulus ignicapilla, Ardea cinerea.

### BC45 | Atlas 1-3

- Richness: 33 → 19 → 10. Flagged interval: 4 gains, 27 losses, 6 persistent species.
- TBI 0.838; null mean 0.665; excess +0.173; Holm p 0.010188.
- Coverage: Good → Average → unavailable.
- Change in Pasture/rangeland: -18.06 percentage points; reference percentile 3.1; local percentile 15.0.
- Starting Annual crops: 0.00%; reference percentile 3.2; local percentile 10.0.
- Starting Pasture/rangeland: 30.84%; reference percentile 96.6; local percentile 100.0.
- Lost species: Anas crecca, Anas platyrhynchos, Charadrius hiaticula, Gallinago gallinago, Pluvialis apricaria, Sternula albifrons, Vanellus vanellus, Acrocephalus schoenobaenus, Acrocephalus scirpaceus, Corvus cornix, Curruca curruca, Cyanistes caeruleus, Emberiza schoeniclus, Erithacus rubecula, Hippolais icterina, Lophophanes cristatus, Motacilla alba, Parus major, Periparus ater, Phylloscopus trochilus, Pica pica, Prunella modularis, Regulus regulus, Saxicola rubetra, Sylvia atricapilla, Turdus philomelos, Athene noctua.
- Gained species: Numenius arquata, Grus grus, Emberiza citrinella, Saxicola rubicola.

### FG85 | Atlas 1-3

- Richness: 54 → 71 → 41. Flagged interval: 13 gains, 26 losses, 28 persistent species.
- TBI 0.582; null mean 0.373; excess +0.209; Holm p 0.010188.
- Coverage: Good → Good → unavailable.
- Change in Broadleaf deciduous forest: 6.95 percentage points; reference percentile 98.4; local percentile 95.0.
- Change in Urban: 4.96 percentage points; reference percentile 92.8; local percentile 90.0.
- Starting Broadleaf deciduous forest: 5.00%; reference percentile 90.2; local percentile 70.0.
- Lost species: Anas platyrhynchos, Acrocephalus scirpaceus, Aegithalos caudatus, Alauda arvensis, Anthus trivialis, Certhia familiaris, Coccothraustes coccothraustes, Curruca communis, Delichon urbicum, Emberiza calandra, Erithacus rubecula, Ficedula hypoleuca, Garrulus glandarius, Hirundo rustica, Linaria cannabina, Loxia curvirostra, Muscicapa striata, Periparus ater, Phylloscopus sibilatrix, Pica pica, Pyrrhula pyrrhula, Regulus regulus, Sturnus vulgaris, Turdus philomelos, Dendrocopos major, Podiceps grisegena.
- Gained species: Buteo buteo, Apus apus, Charadrius dubius, Vanellus vanellus, Streptopelia decaocto, Cuculus canorus, Perdix perdix, Fulica atra, Gallinula chloropus, Acrocephalus palustris, Corvus corax, Corvus frugilegus, Asio otus.
