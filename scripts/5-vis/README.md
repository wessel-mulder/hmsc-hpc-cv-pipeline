# 5 Other

This folder is for release-relevant scripts that do not fit cleanly into preprocessing, HMSC preparation, HMSC fitting, GWR fitting, or figure generation.

It is currently empty. Avoid using this as a general scratch folder. If a script is exploratory, deprecated, or only useful as private provenance, keep it outside the publication-facing workflow or document it clearly before adding it here.

## Inspect the portable S0 dataset

```sh
Rscript data-shared/visualise-s0-species.R
```

This reads only the shared S0 RData and writes a PNG with three atlas presence
maps, the focal species' guild/migration/STI and taxonomic distances to willow
warbler, blackcap and green woodpecker. Default species: Phylloscopus collybita.
Outputs go to `data-shared/`, with supporting CSVs.
Use `--species=Genus_species`, `--rdata=/path/to/file.RData` or
`--output-dir=/path/to/output` to change inputs/outputs; the function
`visualise_s0_species()` is also available by sourcing the script.
Maps use nominal 5-km squares from projected model-site centres and distinguish
non-included cells from absences. They retain each atlas's actual S0 sample;
no additional common-cell filter is applied. Distances are rank-based tree
units, not genetic distances or divergence times. PNG only.
