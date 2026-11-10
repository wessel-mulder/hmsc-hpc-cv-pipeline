Run from the repository root:

```sh
Rscript scripts/5-vis/richness-atlas-visualisation/map-atlas-richness.R
```

Creates four PNGs here: Atlas 1 recorded richness, Atlas 1 chiffchaff and willow warbler occurrence,
and richness across Atlases 1–3 on a shared species-count scale.

Uses `Y` and `Design` directly from the portable S0 RData consumed by S1,
with all 155 retained model species. Richness is the row sum of binary recorded
occurrences, not a fitted model prediction or a count of every raw atlas species.
Chiffchaff is `Phylloscopus_collybita`, and willow warbler is
`Phylloscopus_trochilus`; zero is labelled “Not recorded”.

Retains each atlas's own training cells (1,747 / 1,821 / 1,875), rather than imposing
a new common-cell filter. Therefore, the three-panel figure is descriptive and
does not provide paired cell-level change estimates. The input already applies
the S0 filters: at least five occurrences per species in every atlas, complete
predictors, at least 25% land, and Good/Average Atlas 1–2 coverage (Atlas 3 is
assumed Good). See `scripts/2-hmsc/S0_Data_Definitions.R` for the selection rules.

Uses the temporal-turnover Denmark layout: circular cell-centre points, white
background, grey headers, no axes, and the framed Bornholm inset. Coordinates
are transformed from the bundle's EPSG:23032 to longitude/latitude.
