# A second S2c heatmap: predictor correlations with saved spatial latent factors.
# Correlations are calculated across fitted cells at every retained draw.
s2c_factor_summary <- function(m, meta, interval, targets) {
  n <- length(m$postList[[1]]); nc <- length(m$postList); p <- ncol(m$X)
  variable <- which(apply(m$X, 2, sd) > 0)
  levels <- s2c_spatial_levels(m)
  rows <- list(); dimensions <- list()
  placeholder <- function(level, reason) {
    data.frame(model = meta$model, atlas = meta$atlas, level = level,
      predictor = m$covNames, saved_factor = NA_integer_, n_total = n * nc,
      n_present = 0L, n_finite = 0L, median = NA_real_, lower = NA_real_, upper = NA_real_,
      p_positive = NA_real_, p_negative = NA_real_, interval_excludes_zero = NA,
      rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_, mcse_sd_ratio = NA_real_,
      diagnostic_flag = NA, availability = reason, stringsAsFactors = FALSE)
  }
  if (!length(levels)) return(list(summary = placeholder("", "no spatial factors"), dimensions = data.frame()))
  stopifnot(nrow(m$Pi) == nrow(m$X), !anyNA(m$Pi))
  for (r in levels) {
    level <- s2c_level_name(m, r)
    nf <- vapply(m$postList, function(chain) vapply(chain, function(draw)
      ncol(draw$Eta[[r]]), integer(1)), integer(n))
    max_nf <- max(nf)
    dimensions[[level]] <- data.frame(model = meta$model, atlas = meta$atlas, level = level,
      iteration = rep(seq_len(n), nc), chain = rep(seq_len(nc), each = n), n_factors = as.vector(nf))
    if (!max_nf) {rows[[level]] <- placeholder(level, "no retained spatial factors"); next}
    correlations <- array(NA_real_, c(n, nc, p, max_nf))
    for (k in seq_len(nc)) for (i in seq_len(n)) {
      eta <- m$postList[[k]][[i]]$Eta[[r]]
      stopifnot(all(m$Pi[, r] >= 1), all(m$Pi[, r] <= nrow(eta)), all(is.finite(eta)))
      if (ncol(eta) > 0 && length(variable))
        correlations[i, k, variable, seq_len(ncol(eta))] <-
          suppressWarnings(cor(m$X[, variable, drop = FALSE], eta[m$Pi[, r], , drop = FALSE]))
    }
    factor_rows <- list()
    for (f in seq_len(max_nf)) for (v in seq_len(p)) {
      z <- matrix(correlations[, , v, f], nrow = n, ncol = nc)
      finite <- z[is.finite(z)]
      bounds <- c((1 - interval) / 2, .5, 1 - (1 - interval) / 2)
      q <- if (length(finite)) quantile(finite, bounds, names = FALSE) else rep(NA_real_, 3)
      diag <- s2c_diagnose(z, targets)
      n_present <- sum(nf >= f)
      # Full-chain diagnostics are unavailable for factors absent from some draws.
      bad <- !is.finite(diag$rhat) || diag$rhat >= targets$rhat ||
        !is.finite(diag$ess_bulk) || diag$ess_bulk < targets$ess_per_chain * nc ||
        !is.finite(diag$ess_tail) || diag$ess_tail < targets$ess_per_chain * nc ||
        !is.finite(diag$mcse_sd_ratio) || diag$mcse_sd_ratio > targets$mcse_sd_ratio
      availability <- if (!(v %in% variable)) "constant predictor" else
        if (!length(finite)) "undefined correlation" else
          if (n_present < n * nc) "factor absent in some draws" else "available"
      factor_rows[[length(factor_rows) + 1L]] <- data.frame(model = meta$model, atlas = meta$atlas,
        level = level, predictor = m$covNames[v], saved_factor = f,
        n_total = n * nc, n_present = n_present, n_finite = length(finite),
        median = q[2], lower = q[1], upper = q[3],
        p_positive = if (length(finite)) mean(finite > 0) else NA_real_,
        p_negative = if (length(finite)) mean(finite < 0) else NA_real_,
        interval_excludes_zero = if (length(finite)) q[1] > 0 || q[3] < 0 else NA,
        rhat = diag$rhat, ess_bulk = diag$ess_bulk, ess_tail = diag$ess_tail,
        mcse_sd_ratio = diag$mcse_sd_ratio, diagnostic_flag = if (length(finite)) bad else NA,
        availability = availability, stringsAsFactors = FALSE)
    }
    rows[[level]] <- s2c_bind(factor_rows)
  }
  list(summary = s2c_bind(rows), dimensions = s2c_bind(dimensions))
}

s2c_factor_heatmap <- function(summary, models, interval, atlas_labels, output_file) {
  x <- summary
  predictor_order <- unique(x$predictor)
  nice <- c("(Intercept)" = "Intercept (constant)", tmean_breeding = "Temperature",
    prec_breeding = "Precipitation", hh = "Land-use heterogeneity", perc_urban = "Urban cover",
    perc_cropland = "Cropland cover", perc_pasture = "Pasture cover", perc_grass_shrub = "Grass/shrubland cover",
    perc_forest_deciduous = "Deciduous forest cover", perc_forest_evergreen = "Evergreen forest cover",
    perc_forest_mixed_other = "Mixed/other forest cover")
  predictor_labels <- ifelse(predictor_order %in% names(nice), nice[predictor_order], predictor_order)
  blocks <- unique(x[, c("model", "atlas", "level")])
  blocks <- blocks[order(match(blocks$model, models$model), blocks$level), ]
  block_labels <- vapply(seq_len(nrow(blocks)), function(i) {
    b <- blocks[i, ]
    atlas <- if (b$atlas %in% names(atlas_labels)) atlas_labels[[b$atlas]] else b$atlas
    z <- x[x$model == b$model & x$level == b$level, ]
    label <- if (all(is.na(z$saved_factor))) "No spatial\nfactors" else {
      count <- range(z$n_present / z$n_total)
      paste0("Spatial nf", max(z$saved_factor, na.rm = TRUE), if (min(count) < 1) " (varies)" else "")
    }
    paste(c(label, atlas, if (nzchar(b$level)) b$level), collapse = "\n")
  }, "")
  duplicates <- duplicated(block_labels) | duplicated(block_labels, fromLast = TRUE)
  block_labels[duplicates] <- paste(block_labels[duplicates],
    sub("_Atlas.*", "", blocks$model[duplicates]), sep = "\n")
  block_keys <- paste(blocks$model, blocks$level, sep = " | ")
  x$block <- factor(paste(x$model, x$level, sep = " | "), levels = block_keys, labels = block_labels)
  x$column_id <- paste(match(paste(x$model, x$level, sep = " | "), block_keys),
                       ifelse(is.na(x$saved_factor), 0L, x$saved_factor), sep = "_")
  columns <- unique(x[, c("block", "column_id", "saved_factor", "n_present", "n_total")])
  columns <- columns[order(columns$block, columns$saved_factor, na.last = FALSE), ]
  columns$label <- ifelse(is.na(columns$saved_factor), "None", paste0("LF", columns$saved_factor,
    ifelse(columns$n_present < columns$n_total, sprintf("\n%.0f%% draws", 100 * columns$n_present / columns$n_total), "")))
  x$column_id <- factor(x$column_id, levels = columns$column_id)
  x$predictor_label <- factor(x$predictor, levels = rev(predictor_order), labels = predictor_labels[rev(seq_along(predictor_order))])
  displayed <- round(x$median, 2); displayed[!is.na(displayed) & displayed == 0] <- 0
  x$label <- ifelse(is.na(x$median), ifelse(x$availability == "constant predictor", "const.", "—"),
    paste0(sprintf("%.2f", displayed), ifelse(x$interval_excludes_zero, "*", ""),
           ifelse(x$diagnostic_flag, "†", "")))
  p <- ggplot2::ggplot(x, ggplot2::aes(column_id, predictor_label, fill = median)) +
    ggplot2::geom_tile(colour = "white", linewidth = .6) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, colour = "#202020") +
    ggplot2::facet_grid(~block, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(position = "top", labels = setNames(columns$label, columns$column_id),
                             expand = ggplot2::expansion(add = .02)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = .02)) +
    ggplot2::scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#D98252", midpoint = 0,
      limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1), na.value = "#E5E5E5") +
    ggplot2::labs(title = "Environmental variables and spatial latent factors",
      subtitle = "Posterior median Pearson correlations across all fitted cells; factor columns expand separately for each model.",
      caption = paste0("* ", interval * 100, "% posterior interval excludes zero; † correlation-draw diagnostic targets missed or undefined. Markers are descriptive, not cell-independence significance tests.\n",
        "Factor signs and bases are conditional on the saved fit; LF1 in different models need not describe the same spatial pattern. No new factor alignment is applied.\n",
        "Grey: no spatial factor, constant predictor, or undefined correlation. Draw coverage is labelled for intermittent factors; their intervals condition on presence."),
      x = NULL, y = NULL, fill = "Correlation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
      axis.text.x = ggplot2::element_text(size = 9), axis.text.y = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      panel.spacing.x = grid::unit(.65, "lines"), plot.title = ggplot2::element_text(face = "bold", size = 18),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9), plot.margin = ggplot2::margin(12, 15, 12, 12))
  ggplot2::ggsave(output_file, p, width = max(14, nrow(columns) * .72 + 4),
    height = max(7, length(predictor_order) * .45 + 3), dpi = 200, bg = "white", limitsize = FALSE)
  invisible(x)
}

run_spatial_factor_correlations <- function(output_root, interval = .95,
                                            atlas_labels = c(Atlas1 = "1970s", Atlas2 = "1990s", Atlas3 = "2010s")) {
  stopifnot(length(interval) == 1L, is.finite(interval), interval > 0, interval < 1)
  comparison <- readRDS(file.path(output_root, "comparison.rds"))
  manifest <- comparison$manifest; inputs <- manifest$inputs
  if (!identical(unname(tools::md5sum(names(inputs))), unname(inputs))) stop("Fitted inputs changed since convergence comparison.")
  cache_file <- file.path(output_root, "spatial-factor-correlations.rds")
  result <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
  settings <- list(interval = interval, targets = manifest$targets, schema = 1L)
  reuse <- !is.null(result) && identical(result$inputs, inputs) && identical(result$settings, settings)
  if (!reuse) {
    summaries <- list(); dimensions <- list()
    for (i in seq_len(nrow(manifest$models))) {
      meta <- manifest$models[i, ]; message("Spatial-factor correlations: ", meta$model)
      m <- s2c_read(meta$input)
      s2c_validate(m, meta$samples, meta$thin, meta$n_chains)
      z <- s2c_factor_summary(m, as.list(meta), interval, manifest$targets)
      summaries[[i]] <- z$summary; dimensions[[i]] <- z$dimensions
      rm(m, z); gc(verbose = FALSE)
    }
    result <- list(summary = s2c_bind(summaries), dimensions = s2c_bind(dimensions),
                   inputs = inputs, settings = settings)
  } else message("Reusing validated spatial-factor correlations")
  date_tag <- s2c_date_tag(manifest$models$model)
  file <- file.path(output_root, paste0(date_tag, "-spatial-factor-correlations.png"))
  result$figure_cells <- s2c_factor_heatmap(result$summary, manifest$models, interval, atlas_labels, file)
  if (!identical(unname(tools::md5sum(names(inputs))), unname(inputs))) stop("Fitted inputs changed during spatial correlations.")
  saveRDS(result, cache_file)
  write.csv(result$summary, file.path(output_root, "spatial-factor-correlations.csv"), row.names = FALSE, na = "NA")
  message("Saved spatial-factor comparison: ", file)
  invisible(file)
}
