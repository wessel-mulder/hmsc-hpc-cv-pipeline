# Trait redundancy checks for the model-ready HMSC trait table.
#
# Purpose:
#   Test whether the traits used in the HMSC trait formula are strongly
#   associated with each other before interpreting them as independent
#   ecological axes.
#
# Inputs:
#   Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData
#
# Outputs:
#   notebooks/exploratory/outputs/trait-redundancy/
#     - trait-redundancy-pairwise-tests.csv
#     - trait-counts.csv
#     - sti-by-migration.csv
#     - sti-by-foraging-guild.csv
#     - migration-by-foraging-guild.csv
#     - trait-redundancy-diagnostic.png

rm(list = ls())

input_file <- "Data/preprocessed_data_min_occs_is_5_coverage_is_good_and_average.RData"
output_dir <- "notebooks/exploratory/outputs/trait-redundancy"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the exact preprocessed objects used by S1_Model_Definitions.R.
# The object Tr is a species-by-traits data frame already filtered to the
# modelled species set for min_occs = 5 and Good/Average coverage.
model_data <- new.env()
load(input_file, envir = model_data)

traits <- model_data$Tr
traits$species <- rownames(traits)

# Give short names to the three traits so the analysis code is easy to read.
names(traits)[names(traits) == "Migration_a3_DOF"] <- "migration"
names(traits)[names(traits) == "foraging_guild_consensus"] <- "foraging_guild"
names(traits)[names(traits) == "species_thermal_index"] <- "sti"

traits$migration <- factor(
  traits$migration,
  levels = c(
    "sedentary",
    "sedentary and short-distance",
    "short-distance",
    "short-and long-distance",
    "long-distance"
  )
)
traits$foraging_guild <- factor(traits$foraging_guild)
traits$sti <- as.numeric(traits$sti)

# Keep only complete cases for these three traits. This should retain all
# species in the current model-ready object, but the explicit filter makes the
# script robust if future trait tables contain missing values.
traits_complete <- traits[complete.cases(traits[, c("migration", "foraging_guild", "sti")]), ]

# ---- Small helper functions -------------------------------------------------

format_p <- function(p_value) {
  if (is.na(p_value)) return(NA_character_)
  if (p_value < 0.001) return("<0.001")
  sprintf("%.3f", p_value)
}

eta_squared_oneway <- function(model) {
  aov_table <- anova(model)
  sum_sq <- aov_table[["Sum Sq"]]
  sum_sq[1] / sum(sum_sq)
}

epsilon_squared_kruskal <- function(kruskal_result, n, k) {
  # Epsilon-squared gives a simple non-parametric effect size for the
  # Kruskal-Wallis test. Values near zero indicate little separation among
  # groups; larger values indicate stronger group-level differences.
  as.numeric((kruskal_result$statistic - k + 1) / (n - k))
}

cramers_v <- function(tab) {
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  min_dim <- min(nrow(tab) - 1, ncol(tab) - 1)
  sqrt(as.numeric(test$statistic) / (n * min_dim))
}

# ---- Continuous trait versus categorical traits ----------------------------

sti_migration_lm <- lm(sti ~ migration, data = traits_complete)
sti_migration_kw <- kruskal.test(sti ~ migration, data = traits_complete)

sti_guild_lm <- lm(sti ~ foraging_guild, data = traits_complete)
sti_guild_kw <- kruskal.test(sti ~ foraging_guild, data = traits_complete)

# ---- Categorical trait versus categorical trait ----------------------------

migration_by_guild <- table(traits_complete$migration, traits_complete$foraging_guild)

# There are many foraging guilds, so some expected cell counts are small.
# A simulated chi-square p-value is more appropriate than relying on the
# asymptotic approximation or trying an exact Fisher test on a wide table.
migration_guild_chisq <- suppressWarnings(
  chisq.test(migration_by_guild, simulate.p.value = TRUE, B = 9999)
)

# ---- Summaries --------------------------------------------------------------

pairwise_tests <- data.frame(
  comparison = c(
    "Species thermal index ~ migration",
    "Species thermal index ~ foraging guild",
    "Migration x foraging guild"
  ),
  test = c(
    "Kruskal-Wallis plus one-way linear model effect size",
    "Kruskal-Wallis plus one-way linear model effect size",
    "Chi-square test with simulated p-value"
  ),
  n_species = c(
    nrow(traits_complete),
    nrow(traits_complete),
    sum(migration_by_guild)
  ),
  statistic = c(
    as.numeric(sti_migration_kw$statistic),
    as.numeric(sti_guild_kw$statistic),
    as.numeric(migration_guild_chisq$statistic)
  ),
  df = c(
    as.numeric(sti_migration_kw$parameter),
    as.numeric(sti_guild_kw$parameter),
    as.numeric(migration_guild_chisq$parameter)
  ),
  p_value = c(
    sti_migration_kw$p.value,
    sti_guild_kw$p.value,
    migration_guild_chisq$p.value
  ),
  p_value_formatted = vapply(
    c(sti_migration_kw$p.value, sti_guild_kw$p.value, migration_guild_chisq$p.value),
    format_p,
    character(1)
  ),
  effect_size_type = c(
    "eta_squared_lm / epsilon_squared_kw",
    "eta_squared_lm / epsilon_squared_kw",
    "cramers_v"
  ),
  effect_size_primary = c(
    eta_squared_oneway(sti_migration_lm),
    eta_squared_oneway(sti_guild_lm),
    cramers_v(migration_by_guild)
  ),
  effect_size_secondary = c(
    epsilon_squared_kruskal(
      sti_migration_kw,
      n = nrow(traits_complete),
      k = nlevels(traits_complete$migration)
    ),
    epsilon_squared_kruskal(
      sti_guild_kw,
      n = nrow(traits_complete),
      k = nlevels(traits_complete$foraging_guild)
    ),
    NA_real_
  ),
  adjusted_r2_lm = c(
    summary(sti_migration_lm)$adj.r.squared,
    summary(sti_guild_lm)$adj.r.squared,
    NA_real_
  )
)

trait_counts <- rbind(
  data.frame(
    trait = "migration",
    level = names(sort(table(traits_complete$migration), decreasing = TRUE)),
    n_species = as.integer(sort(table(traits_complete$migration), decreasing = TRUE))
  ),
  data.frame(
    trait = "foraging_guild",
    level = names(sort(table(traits_complete$foraging_guild), decreasing = TRUE)),
    n_species = as.integer(sort(table(traits_complete$foraging_guild), decreasing = TRUE))
  )
)

sti_by_migration <- aggregate(
  sti ~ migration,
  data = traits_complete,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x)
  )
)
sti_by_migration <- do.call(data.frame, sti_by_migration)

sti_by_guild <- aggregate(
  sti ~ foraging_guild,
  data = traits_complete,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x)
  )
)
sti_by_guild <- do.call(data.frame, sti_by_guild)
sti_by_guild <- sti_by_guild[order(-sti_by_guild$sti.n, sti_by_guild$foraging_guild), ]

write.csv(
  pairwise_tests,
  file.path(output_dir, "trait-redundancy-pairwise-tests.csv"),
  row.names = FALSE
)
write.csv(
  trait_counts,
  file.path(output_dir, "trait-counts.csv"),
  row.names = FALSE
)
write.csv(
  sti_by_migration,
  file.path(output_dir, "sti-by-migration.csv"),
  row.names = FALSE
)
write.csv(
  sti_by_guild,
  file.path(output_dir, "sti-by-foraging-guild.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame.matrix(migration_by_guild),
  file.path(output_dir, "migration-by-foraging-guild.csv")
)

# ---- PNG diagnostic ---------------------------------------------------------

# Keep the diagnostic in base graphics so the script does not depend on any
# packages beyond base R. The PNG output follows the project instruction to use
# PNGs for plots unless another format is explicitly requested.
migration_plot_labels <- c(
  "sedentary",
  "sedentary\n+ short",
  "short",
  "short\n+ long",
  "long"
)

png(
  file.path(output_dir, "trait-redundancy-diagnostic.png"),
  width = 3000,
  height = 2200,
  res = 260
)

layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1.1, 1.4))
par(mar = c(6, 5, 3, 1), las = 2)

boxplot(
  sti ~ migration,
  data = traits_complete,
  names = migration_plot_labels,
  col = "grey85",
  border = "grey35",
  ylab = "Species thermal index",
  xlab = "",
  main = "STI by migration strategy",
  cex.axis = 0.8
)

guild_counts <- sort(table(traits_complete$foraging_guild), decreasing = TRUE)
top_guilds <- names(guild_counts)[seq_len(min(12, length(guild_counts)))]
boxplot(
  sti ~ foraging_guild,
  data = traits_complete[traits_complete$foraging_guild %in% top_guilds, ],
  col = "grey85",
  border = "grey35",
  ylab = "Species thermal index",
  xlab = "",
  main = "STI by common foraging guilds",
  cex.axis = 0.65
)

par(mar = c(8, 11, 4, 4), las = 2)
heatmap_matrix <- migration_by_guild[, top_guilds, drop = FALSE]
heatmap_row_labels <- c(
  "sedentary",
  "sedentary + short",
  "short",
  "short + long",
  "long"
)
image(
  x = seq_len(ncol(heatmap_matrix)),
  y = seq_len(nrow(heatmap_matrix)),
  z = t(heatmap_matrix),
  axes = FALSE,
  col = hcl.colors(20, "YlOrRd", rev = FALSE),
  xlab = "",
  ylab = "",
  main = "Migration by common foraging guilds"
)
axis(1, at = seq_len(ncol(heatmap_matrix)), labels = colnames(heatmap_matrix), cex.axis = 0.75)
axis(2, at = seq_len(nrow(heatmap_matrix)), labels = heatmap_row_labels, cex.axis = 0.75)
for (i in seq_len(ncol(heatmap_matrix))) {
  for (j in seq_len(nrow(heatmap_matrix))) {
    text(i, j, labels = heatmap_matrix[j, i], cex = 0.8)
  }
}

dev.off()

print(pairwise_tests)
cat("\nTrait counts written to:", output_dir, "\n")
