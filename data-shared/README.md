# Portable Hmsc inputs for HPC

This folder is tracked in Git. The local `Data/` and `data-new/` folders are not.

- `preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData`:
  S0 objects X, Y, Tr, Design and phy; 5,443 surveys, 155 species and 10 predictors.
  `phy` is the full IOC 15.1 taxonomic surrogate tree (226 tips); S1 prunes it
  to the selected model species. Rank distances are not evolutionary times.

S0 refreshes the shared RData whenever it runs successfully. The source data
and intermediate tables stay in the local `data-new/` directory.

On the HPC, pull this repository, install/load R packages Hmsc and ape,
then run from the repository root:

```sh
git pull
Rscript scripts/2-hmsc/S1_Model_Definitions.R
```

S1 reads only the RData in this folder, including its bundled tree. Its current settings create
separate unfitted models for Atlases 1, 2 and 3. It writes model definitions,
script provenance under ignored `HmscOutputs/`. Cross-atlas test files are not generated.
Set `CH1_MODEL_OUTPUT_DIR` to choose another output location. S1 constructs
models; it does not sample or submit an HPC fitting job. It passes raw predictors
to Hmsc with `XScale=TRUE`, storing the training means and standard deviations.
The shared run configuration is `scripts/2-hmsc/run-config.tsv`, with run date
`2026-11-02`. Downstream summaries and predictions use all retained draws by
default. See `scripts/2-hmsc/README.md` for the full HPC command sequence.
Production fitting has not yet been run for this dataset.

The predictors are breeding temperature, breeding precipitation, heterogeneity,
class-22 cropland, urban, pasture, grass/shrub, deciduous forest (43+44), evergreen
forest (42+41), and mixed/other forest (40+45). Cover values are percentages.
Species have >=5 occurrences in every atlas; Good/Average coverage is used for
Atlases 1-2, with the existing assumed-Good Atlas 3 rule. Cells with >75%
ocean + freshwater are excluded. The Corvus cornix STI uses the explicit
Corvus corone proxy. Raw data access is not required merely to run S1.

Local validation: S1 constructed all three atlas model objects successfully;
each retained 155 species and 155 matching tree tips from the full 226-tip tree.

## Visualise the shared input

```sh
Rscript data-shared/visualise-s0-species.R
```

Run from the repository root. The script, PNG summary and supporting CSVs are
stored here beside the S0 RData.
