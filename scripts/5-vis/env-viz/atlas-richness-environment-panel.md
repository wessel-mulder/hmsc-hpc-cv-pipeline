# Atlas richness and environment panel
Run: `Rscript scripts/5-vis/env-viz/atlas-richness-environment-panel.R`
Training input: data-shared/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData
Land-cover input: data-new/0-pre/derived/environmental/X_environment.csv
Same 1697 cells in each atlas; recorded richness of the same 155 model species.
Five columns correspond to the five requested variables; three rows represent the atlas periods.
Temperature and precipitation summarize March–June in original units. Each continuous variable has one scale across all three atlases.
Heterogeneity is mean annual Simpson diversity across all 15 original HILDA classes, including ocean and freshwater.
Dominant land use is the largest period-mean terrestrial cover group; ocean and freshwater excluded. Forest combines classes 40–45.
Groups: cropland (22), forest (40–45), pasture (33), urban (11), grass/shrub (55), tree crops (23), agroforestry (24), sparse/bare (66).
Exact ties use the first group in that order. Dominance need not exceed 50% of a cell. Legend shows groups dominant in at least one cell.
Climate years: 1971–1974, 1993–1996, 2014–2016; HILDA years: 1971–1974, 1993–1996, 2014–2017.
   
    Cropland Forest Pasture Urban Grass / shrub Tree crops Agroforestry
  1     1317    199      66    60            53          0            0
  2     1255    220      66    68            86          0            0
  3     1261    235      51    79            66          0            0
   
    Sparse / bare
  1             2
  2             2
  3             5
