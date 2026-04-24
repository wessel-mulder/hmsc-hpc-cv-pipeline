read_effort_coverage <- function(dir) {
  atlas1 <- readxl::read_excel(
    file.path(dir, "effort", "dækning atlas1+2.xlsx"),
    sheet = "kvadrater atlas1"
  )
  atlas2 <- readxl::read_excel(
    file.path(dir, "effort", "dækning atlas1+2.xlsx"),
    sheet = "kvadrater atlas2"
  )

  atlas2 <- atlas2 |>
    dplyr::select(-UNDS_KODE_) |>
    dplyr::rename(UNDS_KODE_ = UNDS_KOD_1)

  key <- readxl::read_excel(
    file.path(dir, "effort", "dækning atlas1+2.xlsx"),
    sheet = "dækning"
  ) |>
    dplyr::mutate(
      coverage = dplyr::case_when(
        UNDSKODE == 0 ~ "Not understood",
        UNDSKODE == 1 ~ "Good",
        UNDSKODE == 2 ~ "Average",
        UNDSKODE == 3 ~ "Bad"
      )
    )

  effort <- rbind(
    atlas1 |> dplyr::mutate(atlas = 1),
    atlas2 |> dplyr::mutate(atlas = 2)
  ) |>
    dplyr::left_join(
      key |> dplyr::rename(UNDS_KODE_ = UNDSKODE),
      by = "UNDS_KODE_"
    ) |>
    dplyr::rename(site = KVADRATKOD)

  effort$coverage <- factor(
    effort$coverage,
    levels = c("Good", "Average", "Bad", "Not understood")
  )

  effort |>
    dplyr::mutate(survey = paste(site, atlas, sep = "_")) |>
    dplyr::select(survey, coverage)
}

lulc_lookup <- c(
  "0" = "ocean",
  "11" = "urban",
  "22" = "cropland",
  "33" = "pasture",
  "44" = "forest",
  "55" = "grass_shrub",
  "66" = "other",
  "77" = "water"
)

clean_lulc_columns <- function(df) {
  df |>
    dplyr::rename_with(
      .cols = tidyselect::matches("^LULC_"),
      .fn = function(col) {
        lulc_code <- stringr::str_extract(col, "\\d+")
        paste0("perc_", lulc_lookup[lulc_code])
      }
    ) |>
    dplyr::mutate(
      perc_fresh_saltwater = perc_water + perc_ocean
    ) |>
    dplyr::select(-perc_other, -perc_ocean, -perc_water)
}

sti_taxonomy_map <- c(
  "Accipiter_gentilis" = "Astur_gentilis",
  "Charadrius_alexandrinus" = "Anarhynchus_alexandrinus",
  "Larus_ridibundus" = "Chroicocephalus_ridibundus",
  "Corvus_monedula" = "Coloeus_monedula",
  "Sylvia_communis" = "Curruca_communis",
  "Sylvia_curruca" = "Curruca_curruca"
)

standardise_sti_species_names <- function(sti) {
  sti$SPECIES <- gsub(" ", "_", sti$SPECIES)
  sti$SPECIES <- ifelse(
    sti$SPECIES %in% names(sti_taxonomy_map),
    sti_taxonomy_map[sti$SPECIES],
    sti$SPECIES
  )

  if ("Corvus_corone" %in% sti$SPECIES && !"Corvus_cornix" %in% sti$SPECIES) {
    corone_data <- sti[sti$SPECIES == "Corvus_corone", ]
    corone_data$SPECIES <- "Corvus_cornix"
    sti <- rbind(sti, corone_data)
  }

  sti
}
