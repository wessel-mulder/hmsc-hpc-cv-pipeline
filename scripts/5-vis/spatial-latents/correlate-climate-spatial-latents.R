# Run from the repository root:
# Rscript scripts/5-vis/spatial-latents/correlate-climate-spatial-latents.R
# Descriptive Pearson correlations across cells, conditional on fitted models.
# No independence-based correlation tests or p-values are calculated.

required <- c("ggplot2", "posterior")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))
model_stamp <- "2026-11-02_06-55-21"
model_root <- file.path("data-new", "2-hmsc", "derived", "models")
output_dir <- file.path("scripts", "5-vis", "spatial-latents")
climate_names <- c("tmean_breeding", "prec_breeding")
climate_labels <- c("Temperature", "Precipitation")
atlas_labels <- c("1970s", "1990s", "2010s")

read_model <- function(atlas) {
  folder <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas)
  path <- file.path(model_root, folder, "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata")
  e <- new.env(parent = emptyenv()); load(path, envir = e)
  m <- e$fitted_model$posteriors
  if (length(m$postList) != 4L || any(lengths(m$postList) != 250L) || m$nr != 1L) {
    stop("Unexpected model dimensions: ", folder)
  }
  sites <- as.character(m$studyDesign$site)
  if (anyDuplicated(sites) || !identical(rownames(m$XData), rownames(m$studyDesign)) ||
      !identical(as.character(m$dfPi$site), sites) ||
      !identical(as.integer(m$Pi[, 1]), match(sites, levels(m$dfPi$site)))) {
    stop("Site order mismatch: ", folder)
  }
  climate <- as.matrix(m$XData[, climate_names, drop = FALSE])
  if (any(!is.finite(climate)) || any(apply(climate, 2, sd) <= 0)) stop("Invalid climate values.")
  eta <- array(NA_real_, c(length(sites), 5L, 250L, 4L))
  for (chain in 1:4) for (iteration in 1:250) {
    draw <- m$postList[[chain]][[iteration]]$Eta[[1]]
    if (!identical(dim(draw), c(nlevels(m$dfPi$site), 5L)) || any(!is.finite(draw))) {
      stop("Invalid Eta draw.")
    }
    eta[, , iteration, chain] <- draw[m$Pi[, 1], , drop = FALSE]
  }
  list(sites = sites, climate = climate, eta = eta)
}
models <- lapply(1:3, read_model)
stopifnot(identical(vapply(models, function(m) length(m$sites), integer(1)), c(1747L, 1821L, 1875L)))
common_sites <- sort(Reduce(intersect, lapply(models, `[[`, "sites")))
stopifnot(length(common_sites) == 1697L)

summary_rows <- list(); draw_rows <- list(); chain_rows <- list()
row_i <- 0L
for (atlas in 1:3) for (cell_set in c("fitted", "common")) {
  m <- models[[atlas]]
  ix <- if (cell_set == "fitted") seq_along(m$sites) else match(common_sites, m$sites)
  stopifnot(!anyNA(ix))
  x <- m$climate[ix, , drop = FALSE]
  mean_eta <- apply(m$eta[ix, , , , drop = FALSE], c(1, 2), mean)
  mean_map_r <- cor(x, mean_eta, method = "pearson")
  r <- array(NA_real_, c(250L, 4L, 2L, 5L))
  for (chain in 1:4) for (iteration in 1:250) {
    r[iteration, chain, , ] <- cor(x, m$eta[ix, , iteration, chain], method = "pearson")
  }
  if (any(!is.finite(r)) || any(abs(r) > 1 + 1e-12)) stop("Invalid correlations.")
  # Verify vectorised correlations independently for one draw, both climates/all factors.
  for (cl in 1:2) for (fac in 1:5) {
    direct <- cor(x[, cl], m$eta[ix, fac, 1, 1])
    stopifnot(abs(direct - r[1, 1, cl, fac]) < 1e-12)
    z <- r[, , cl, fac]
    q <- quantile(z, c(0.025, 0.5, 0.975), names = FALSE)
    row_i <- row_i + 1L
    key <- data.frame(atlas = atlas, cell_set = cell_set, n_cells = length(ix),
                      climate = climate_labels[[cl]], factor = fac)
    summary_rows[[row_i]] <- cbind(key, data.frame(
      mean_map_r = mean_map_r[cl, fac], median_r = q[[2]], lower_r = q[[1]], upper_r = q[[3]],
      pr_positive = mean(z > 0), rhat = posterior::rhat(z),
      ess_bulk = posterior::ess_bulk(z), ess_tail = posterior::ess_tail(z),
      mcse_mean = posterior::mcse_mean(z)))
    draw_rows[[row_i]] <- cbind(key[rep(1L, 1000L), ], data.frame(
      chain = rep(1:4, each = 250L), iteration = rep(1:250, times = 4L), r = as.vector(z)))
    chain_rows[[row_i]] <- cbind(key[rep(1L, 4L), ], data.frame(
      chain = 1:4, mean_r = colMeans(z), median_r = apply(z, 2, median)))
  }
}
correlation_summary <- do.call(rbind, summary_rows)
correlation_draws <- do.call(rbind, draw_rows)
correlation_chains <- do.call(rbind, chain_rows)
rownames(correlation_summary) <- rownames(correlation_draws) <- rownames(correlation_chains) <- NULL
stopifnot(nrow(correlation_summary) == 60L, nrow(correlation_draws) == 60000L)
if (any(!is.finite(as.matrix(correlation_summary[, c("rhat", "ess_bulk", "ess_tail", "mcse_mean")])))) {
  stop("Non-finite correlation diagnostics; inspect draws before summarising.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (cell_set in c("fitted", "common")) {
  dat <- correlation_summary[correlation_summary$cell_set == cell_set, ]
  dat$atlas <- factor(dat$atlas, levels = 1:3, labels = paste0("Atlas ", 1:3, " (", atlas_labels, ")"))
  dat$climate <- factor(dat$climate, levels = climate_labels)
  dat$factor <- factor(dat$factor, levels = 5:1, labels = paste("Factor", 5:1))
  points <- rbind(data.frame(dat, estimate = dat$median_r, summary_type = "Posterior median"),
                  data.frame(dat, estimate = dat$mean_map_r, summary_type = "Mean-map correlation"))
  p <- ggplot(dat, aes(y = factor)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65") +
    geom_segment(aes(x = lower_r, xend = upper_r, yend = factor), linewidth = 0.7, colour = "#2C7FB8") +
    geom_point(data = points, aes(x = estimate, shape = summary_type), size = 2.5, colour = "#2C7FB8") +
    facet_grid(atlas ~ climate) +
    scale_shape_manual(values = c("Posterior median" = 16, "Mean-map correlation" = 5)) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    labs(title = "Climate correlations with spatial latent factors",
         subtitle = if (cell_set == "fitted") "All fitted cells; 95% pooled-draw intervals (chain disagreement limits inference)"
                    else "Same 1,697 cells; 95% pooled-draw intervals (chain disagreement limits inference)",
         x = "Pearson correlation across cells", y = NULL, shape = NULL,
         caption = "All 1,000 draws; saved within-atlas factor alignment. Mean-map correlation differs from the draw-level median.\nMany derived correlations have poor R-hat: inspect the report diagnostics before interpreting intervals. No cell-independence p-values.") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
  path <- file.path(output_dir, paste0("climate-spatial-latent-correlations-", cell_set, ".png"))
  ggsave(path, p, width = 12, height = 10, dpi = 300, bg = "white")
  message("Saved: ", path)
}
chain_dat <- correlation_chains[correlation_chains$cell_set == "fitted", ]
chain_dat$atlas <- factor(chain_dat$atlas, levels = 1:3, labels = paste0("Atlas ", 1:3, " (", atlas_labels, ")"))
chain_dat$climate <- factor(chain_dat$climate, levels = climate_labels)
chain_dat$chain <- factor(chain_dat$chain)
chain_plot <- ggplot(chain_dat, aes(factor, mean_r, colour = chain, group = chain)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65") +
  geom_line(linewidth = 0.7) + geom_point(size = 2) +
  facet_grid(atlas ~ climate) +
  scale_x_continuous(breaks = 1:5) + scale_y_continuous(limits = c(-1, 1)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  labs(title = "Climate–latent correlations: agreement among chains",
       subtitle = "Mean draw-level correlation within each chain; all fitted cells",
       x = "Within-atlas factor", y = "Chain mean Pearson correlation", colour = "Chain",
       caption = "Lines join factor labels for readability; factor order is not a continuous ecological gradient.\nDifferent chain means can reflect incomplete mixing and/or unresolved latent-factor alignment.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
ggsave(file.path(output_dir, "climate-spatial-latent-correlations-chains.png"), chain_plot,
       width = 12, height = 10, dpi = 300, bg = "white")
message("Validated 60 climate-factor summaries, including 1,697-cell sensitivity; 60,000 draw-level correlations.")
message("Fitted-cell correlations with R-hat >= 1.1: ",
        sum(correlation_summary$rhat[correlation_summary$cell_set == "fitted"] >= 1.1), "/30.")
