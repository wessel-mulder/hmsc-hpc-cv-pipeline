# S2c helpers only. No dependency on or changes to S3's global settings.
s2c_bind <- function(xs) {
  xs <- Filter(function(x) !is.null(x) && nrow(x) > 0, xs)
  if (!length(xs)) return(data.frame())
  keys <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(x) {
    for (v in setdiff(keys, names(x))) x[[v]] <- NA
    x[, keys, drop = FALSE]
  }))
}
s2c_stat <- function(fun, z, ...) {
  tryCatch(suppressWarnings(as.numeric(fun(z, ...))), error = function(e) NA_real_)
}
s2c_diagnose <- function(z, targets) {
  stopifnot(is.matrix(z), nrow(z) > 0, ncol(z) > 0)
  finite <- all(is.finite(z))
  usable <- finite && nrow(z) >= 4 && ncol(z) >= 2
  s <- if (finite) sd(as.vector(z)) else NA_real_
  q <- if (finite) quantile(z, c(.025, .5, .975), names = FALSE) else rep(NA_real_, 3)
  d <- if (usable) c(s2c_stat(posterior::rhat, z),
                      s2c_stat(posterior::ess_bulk, z),
                      s2c_stat(posterior::ess_tail, z),
                      s2c_stat(posterior::mcse_mean, z)) else rep(NA_real_, 4)
  endpoint <- if (usable) s2c_stat(posterior::mcse_quantile, z, probs = c(.025, .975)) else NA_real_
  if (length(endpoint) != 2) endpoint <- rep(NA_real_, 2)
  ratio <- if (is.finite(s) && s > 0) d[4] / s else NA_real_
  reason <- if (!finite) "nonfinite draws" else if (!usable) "too few chains/draws" else
    if (s == 0) "constant sampled draws" else if (any(!is.finite(c(d, ratio, endpoint))))
      "undefined diagnostic" else ""
  data.frame(n_draws = length(z), n_chains = ncol(z), median = q[2], sd = s,
    lower = q[1], upper = q[3], rhat = d[1], ess_bulk = d[2], ess_tail = d[3],
    ess_bulk_relative = d[2] / length(z), ess_tail_relative = d[3] / length(z),
    mcse_mean = d[4], mcse_sd_ratio = ratio, mcse_lower = endpoint[1], mcse_upper = endpoint[2],
    unavailable_reason = reason, row.names = NULL)
}
s2c_metrics <- c("rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio")
s2c_target <- function(metric, n_chains, targets) {
  switch(metric, rhat = targets$rhat, ess_bulk = targets$ess_per_chain * n_chains,
         ess_tail = targets$ess_per_chain * n_chains, mcse_sd_ratio = targets$mcse_sd_ratio)
}
s2c_departure <- function(metric, value, target) {
  if (metric %in% c("ess_bulk", "ess_tail")) target - value else value - target
}
s2c_failure <- function(metric, value, target) {
  if (metric == "rhat") value >= target else s2c_departure(metric, value, target) > 0
}
s2c_array <- function(m, getter, size) {
  counts <- lengths(m$postList)
  stopifnot(length(counts) > 0, all(counts == counts[1]))
  a <- array(NA_real_, c(counts[1], length(counts), size))
  for (k in seq_along(m$postList)) for (i in seq_len(counts[1])) {
    z <- getter(m$postList[[k]][[i]])
    if (length(z) != size) stop("Parameter dimensions change across draws.")
    a[i, k, ] <- as.vector(z)
  }
  a
}
s2c_matrix <- function(a, j) matrix(a[, , j], nrow = dim(a)[1], ncol = dim(a)[2])
s2c_key <- function(family, parameter, predictor = "", species = "", level = "") {
  data.frame(family = family, parameter = parameter, predictor = predictor,
             species = species, level = level, stringsAsFactors = FALSE)
}
s2c_summarise <- function(a, key, targets, chains = FALSE, half = FALSE) {
  stopifnot(dim(a)[3] == nrow(key), !anyDuplicated(key$parameter))
  ds <- vector("list", nrow(key)); cs <- vector("list", nrow(key)); hs <- vector("list", nrow(key))
  for (j in seq_len(nrow(key))) {
    z <- s2c_matrix(a, j)
    ds[[j]] <- cbind(key[j, , drop = FALSE], s2c_diagnose(z, targets))
    if (chains) cs[[j]] <- s2c_bind(lapply(seq_len(ncol(z)), function(k) {
      v <- z[, k]; q <- if (all(is.finite(v))) quantile(v, c(.025, .5, .975), names = FALSE) else rep(NA_real_, 3)
      cbind(key[j, , drop = FALSE], data.frame(chain = k, mean = mean(v), sd = sd(v),
                                            lower = q[1], median = q[2], upper = q[3]))
    }))
    if (half) {
      h <- s2c_diagnose(z[seq.int(floor(nrow(z) / 2) + 1L, nrow(z)), , drop = FALSE], targets)
      h$full_median <- ds[[j]]$median; h$full_rhat <- ds[[j]]$rhat
      h$median_shift_in_full_sd <- if (is.finite(ds[[j]]$sd) && ds[[j]]$sd > 0)
        (h$median - ds[[j]]$median) / ds[[j]]$sd else NA_real_
      hs[[j]] <- cbind(key[j, , drop = FALSE], h)
    }
  }
  list(diagnostics = s2c_bind(ds), chains = s2c_bind(cs), half = s2c_bind(hs))
}
s2c_read <- function(input) {
  e <- new.env(parent = emptyenv()); load(input, e)
  m <- e$fitted_model$posteriors
  if (is.null(m$postList) || !length(m$postList)) stop("Missing posterior chains: ", input)
  m
}
s2c_validate <- function(m, samples, thin, chains) {
  stopifnot(length(m$postList) == chains, all(lengths(m$postList) == samples),
            m$samples == samples, m$thin == thin, length(m$covNames) > 0,
            length(m$spNames) > 0, !anyDuplicated(m$covNames), !anyDuplicated(m$spNames),
            is.matrix(m$X), ncol(m$X) == length(m$covNames), all(is.finite(m$X)),
            identical(colnames(m$X), m$covNames), ncol(m$Y) == length(m$spNames),
            identical(colnames(m$Y), m$spNames), nrow(m$X) == nrow(m$Y),
            length(m$studyDesign$site) == nrow(m$X), !anyNA(m$studyDesign$site),
            !anyDuplicated(m$studyDesign$site))
  for (ch in m$postList) for (d in ch) {
    if (!identical(dim(d$Beta), c(length(m$covNames), length(m$spNames))))
      stop("Saved Beta dimensions do not match predictor/species names.")
    if (length(d$Gamma) && !all(is.na(d$Gamma)) &&
        !identical(dim(d$Gamma), c(length(m$covNames), length(m$trNames))))
      stop("Saved Gamma dimensions do not match predictor/trait names.")
    if (length(d$V) && !all(is.na(d$V)) &&
        !identical(dim(d$V), rep(length(m$covNames), 2)))
      stop("Saved V dimensions do not match predictor names.")
    if (length(m$ranLevels) && (length(d$Eta) != length(m$ranLevels) ||
                               length(d$Lambda) != length(m$ranLevels)))
      stop("Saved random levels do not match Eta/Lambda lists.")
  }
}
s2c_spatial_levels <- function(m) which(vapply(m$ranLevels, function(r) isTRUE(r$sDim > 0), logical(1)))
s2c_level_name <- function(m, r) {
  n <- names(m$ranLevels)[r]
  if (!length(n) || is.na(n) || !nzchar(n)) paste0("level", r) else n
}
s2c_date_tag <- function(models) {
  dates <- sort(unique(regmatches(models, regexpr("^[0-9]{4}-[0-9]{2}-[0-9]{2}", models))))
  if (!length(dates)) stop("Selected model folders must start with run dates.")
  if (length(unique(substr(dates, 1, 7))) == 1L)
    paste(c(dates[1], substr(dates[-1], 9, 10)), collapse = "-") else paste(dates, collapse = "_")
}
s2c_metadata <- function(m, input) {
  model <- basename(dirname(dirname(dirname(input))))
  if (!grepl("Atlas[0-9]+", model)) stop("Cannot identify atlas from folder: ", model)
  atlas <- regmatches(model, regexpr("Atlas[0-9]+", model))
  nf <- if (length(m$ranLevels)) vapply(m$postList, function(ch)
    max(vapply(ch, function(d) sum(vapply(d$Lambda, nrow, integer(1))), integer(1))), integer(1)) else 0L
  spatial <- length(s2c_spatial_levels(m)) > 0
  variant <- if (spatial) paste0("Space nf", min(nf), if (max(nf) != min(nf)) paste0("–", max(nf)) else "") else
    if (length(m$ranLevels)) "Nonspatial random effects" else "Environment only"
  # The full folder remains the unique ID; the short label also distinguishes repeated dates/times.
  stamp <- sub("_Atlas.*", "", model)
  list(model = model, atlas = atlas, variant = variant,
       label = paste(variant, stamp, sep = "\n"), factor_order = max(nf), input = normalizePath(input),
       sites = as.character(m$studyDesign$site), species = m$spNames, predictors = m$covNames,
       X = m$X, Y = m$Y, XScalePar = m$XScalePar,
       sampling = c(samples = m$samples, thin = m$thin, chains = length(m$postList)),
       row = data.frame(model = model, atlas = atlas, variant = variant,
         label = paste(variant, stamp, sep = "\n"), factor_order = max(nf),
         n_cells = nrow(m$X), n_species = length(m$spNames), n_predictors = length(m$covNames),
         n_chains = length(m$postList), samples = m$samples, thin = m$thin,
         spatial = spatial, input = normalizePath(input), stringsAsFactors = FALSE))
}
s2c_tag <- function(x, meta) {
  if (!nrow(x)) return(x)
  cbind(meta$row[rep(1L, nrow(x)), c("model", "atlas", "label", "variant"), drop = FALSE], x, row.names = NULL)
}
s2c_parameter_families <- function(m, targets) {
  p <- length(m$covNames); s <- length(m$spNames)
  bk <- expand.grid(predictor = m$covNames, species = m$spNames, stringsAsFactors = FALSE)
  key <- s2c_key("Beta", paste(bk$predictor, bk$species, sep = " | "), bk$predictor, bk$species)
  beta <- s2c_array(m, function(d) d$Beta, p * s)
  out <- s2c_summarise(beta, key, targets, chains = TRUE, half = TRUE)
  add <- function(a, k) {
    z <- s2c_summarise(a, k, targets)
    out$diagnostics <<- s2c_bind(list(out$diagnostics, z$diagnostics))
  }
  first <- m$postList[[1]][[1]]
  if (length(first$Gamma) && !all(is.na(first$Gamma))) {
    k <- expand.grid(predictor = m$covNames, trait = m$trNames, stringsAsFactors = FALSE)
    stopifnot(length(first$Gamma) == nrow(k))
    add(s2c_array(m, function(d) d$Gamma, nrow(k)),
        s2c_key("Gamma", paste(k$predictor, k$trait, sep = " | "), k$predictor))
  }
  if (length(first$V) && !all(is.na(first$V))) {
    ij <- which(upper.tri(first$V, diag = TRUE), arr.ind = TRUE)
    add(s2c_array(m, function(d) d$V[ij], nrow(ij)),
        s2c_key("V", paste(m$covNames[ij[, 1]], m$covNames[ij[, 2]], sep = " | ")))
  }
  if (length(first$rho) && !all(is.na(first$rho))) {
    # A one-point prior support is structural; constant sampled estimates are not.
    fixed <- !is.null(m$rhopw) && sum(m$rhopw[, 2] > 0) == 1L
    if (!fixed) add(s2c_array(m, function(d) d$rho, length(first$rho)),
                   s2c_key("rho", paste0("rho[", seq_along(first$rho), "]")))
  }
  # Probit residual variances are fixed by the response distribution and excluded.
  sigma_index <- if (is.matrix(m$distr)) which(m$distr[, "variance"] != 0) else integer()
  if (length(sigma_index)) add(s2c_array(m, function(d) d$sigma[sigma_index], length(sigma_index)),
                              s2c_key("sigma", m$spNames[sigma_index], species = m$spNames[sigma_index]))
  out$beta <- beta; out$beta_key <- key
  out
}
s2c_spatial_cor <- function(x, eta, lambda) {
  xc <- scale(x); ec <- scale(eta, center = TRUE, scale = FALSE)
  variance <- colSums(lambda * (crossprod(ec) %*% lambda)) / (nrow(x) - 1)
  result <- sweep(crossprod(xc, ec) %*% lambda / (nrow(x) - 1), 2, sqrt(pmax(variance, 0)), "/")
  result[, variance <= 0] <- NA_real_
  result
}
s2c_project <- function(x, eta, lambda) {
  if (qr(x)$rank < ncol(x)) return(NULL)
  qr.coef(qr(x), eta) %*% lambda
}
s2c_spatial <- function(m, beta, beta_key, targets) {
  out <- list(diagnostics = data.frame(), fields = data.frame(), dimensions = data.frame(),
              compensation = data.frame(), notes = character())
  levels <- s2c_spatial_levels(m)
  if (!length(levels)) {out$notes <- "Spatial diagnostics inapplicable: no spatial random levels."; return(out)}
  n <- dim(beta)[1]; nc <- dim(beta)[2]; p <- length(m$covNames); s <- length(m$spNames)
  stopifnot(nrow(m$Pi) == nrow(m$X))
  x <- m$X; variable <- which(apply(x, 2, sd) > 0); fullrank <- qr(x)$rank == p
  for (r in levels) {
    level <- s2c_level_name(m, r)
    message("  Spatial level ", level, ": associations, overlap and compensation")
    ij <- which(upper.tri(matrix(FALSE, s, s)), arr.ind = TRUE)
    # All off-diagonal species pairs, no sampling; process diagnostics in small blocks.
    pairs <- matrix(NA_real_, n * nc, nrow(ij))
    overlap <- array(NA_real_, c(n, nc, length(variable) * s))
    projected <- if (fullrank) array(NA_real_, dim(beta)) else NULL
    fields <- vector("list", n * nc); dimensions <- vector("list", n * nc)
    alpha_list <- vector("list", n * nc)
    for (k in seq_len(nc)) for (i in seq_len(n)) {
      j <- (k - 1L) * n + i; d <- m$postList[[k]][[i]]
      lam <- d$Lambda[[r]]; raw <- d$Eta[[r]]
      stopifnot(ncol(lam) == s, ncol(raw) == nrow(lam), all(is.finite(raw)), all(is.finite(lam)),
                all(m$Pi[, r] >= 1), all(m$Pi[, r] <= nrow(raw)))
      eta <- raw[m$Pi[, r], , drop = FALSE]
      covariance <- crossprod(lam); variance <- diag(covariance)
      correlation <- covariance / sqrt(outer(variance, variance))
      pairs[j, ] <- correlation[ij]
      if (length(variable) && ncol(eta) > 0) {
        overlap[i, k, ] <- as.vector(s2c_spatial_cor(x[, variable, drop = FALSE], eta, lam))
        cr <- suppressWarnings(cor(x[, variable, drop = FALSE], eta))
        fields[[j]] <- data.frame(chain = k, iteration = i,
          predictor = rep(m$covNames[variable], ncol(eta)),
          saved_factor = rep(seq_len(ncol(eta)), each = length(variable)), correlation = as.vector(cr))
      }
      dimensions[[j]] <- data.frame(chain = k, iteration = i, n_factors = ncol(eta))
      if (fullrank) projected[i, k, ] <- as.vector(s2c_project(x, eta, lam))
      alpha <- d$Alpha[[r]]
      if (length(alpha)) {
        grid <- m$ranLevels[[r]]$alphapw
        stopifnot(all(is.finite(alpha)), all(alpha == as.integer(alpha)), all(alpha >= 1), all(alpha <= nrow(grid)))
        alpha_list[[j]] <- grid[alpha, 1] # Hmsc stores grid indices, not distances.
      }
    }
    stable <- length(unique(vapply(alpha_list, length, integer(1)))) == 1 && length(alpha_list[[1]]) > 0
    alpha_fixed <- sum(m$ranLevels[[r]]$alphapw[, 2] > 0) == 1L
    if (stable && !alpha_fixed) {
      a <- array(unlist(alpha_list), c(length(alpha_list[[1]]), n, nc))
      a <- aperm(a, c(2, 3, 1))
      key <- s2c_key("Alpha (saved basis)", paste0(level, " | saved factor ", seq_len(dim(a)[3])), level = level)
      out$diagnostics <- s2c_bind(list(out$diagnostics, s2c_summarise(a, key, targets)$diagnostics))
    } else out$notes <- c(out$notes, paste(level, if (alpha_fixed)
      "Alpha structurally fixed: excluded." else "Alpha diagnostics unavailable: changing/absent factor dimensions."))
    blocks <- split(seq_len(nrow(ij)), ceiling(seq_len(nrow(ij)) / 500))
    for (bb in seq_along(blocks)) {
      sel <- blocks[[bb]]
      a <- array(pairs[, sel, drop = FALSE], c(n, nc, length(sel)))
      key <- s2c_key("Residual correlation", paste(level, m$spNames[ij[sel, 1]], m$spNames[ij[sel, 2]], sep = " | "), level = level)
      out$diagnostics <- s2c_bind(list(out$diagnostics, s2c_summarise(a, key, targets)$diagnostics))
      if (bb %% 8 == 0) message("    Residual associations: ", max(sel), "/", nrow(ij))
    }
    if (length(variable)) {
      bk <- expand.grid(predictor = m$covNames[variable], species = m$spNames, stringsAsFactors = FALSE)
      key <- s2c_key("Spatial overlap", paste(level, bk$predictor, bk$species, sep = " | "), bk$predictor, bk$species, level)
      out$diagnostics <- s2c_bind(list(out$diagnostics, s2c_summarise(overlap, key, targets)$diagnostics))
    }
    f <- s2c_bind(fields)
    if (nrow(f)) {
      groups <- split(f, interaction(f$chain, f$predictor, f$saved_factor, drop = TRUE))
      f <- s2c_bind(lapply(groups, function(z) {
        values <- z$correlation[is.finite(z$correlation)]
        q <- if (length(values)) quantile(values, c(.025, .5, .975), names = FALSE) else rep(NA_real_, 3)
        data.frame(chain = z$chain[1], predictor = z$predictor[1], saved_factor = z$saved_factor[1],
                   n_finite = length(values), n_draws = nrow(z), lower = q[1], median = q[2], upper = q[3])
      }))
      f$level <- level; out$fields <- s2c_bind(list(out$fields, f))
    }
    dims <- s2c_bind(dimensions); dims$level <- level
    out$dimensions <- s2c_bind(list(out$dimensions, dims))
    if (fullrank) {
      for (component in c("Projected space", "Beta + projected space")) {
        a <- if (component == "Projected space") projected else projected + beta
        key <- beta_key; key$family <- component; key$level <- level
        key$parameter <- paste(level, key$parameter, sep = " | ")
        out$diagnostics <- s2c_bind(list(out$diagnostics, s2c_summarise(a, key, targets)$diagnostics))
      }
      comp <- s2c_bind(lapply(seq_len(nrow(beta_key)), function(j) {
        b <- s2c_matrix(beta, j); a <- s2c_matrix(projected, j)
        s2c_bind(lapply(seq_len(nc), function(k) cbind(beta_key[j, , drop = FALSE],
          data.frame(chain = k, correlation = suppressWarnings(cor(b[, k], a[, k]))))))
      }))
      comp$level <- level; out$compensation <- s2c_bind(list(out$compensation, comp))
    } else out$notes <- c(out$notes, paste(level, "Rank-deficient X: projections and compensation skipped."))
    out$notes <- c(out$notes, paste(level, "uses all", nrow(x), "fitted cells; factor signs/basis are conditional on saved draws."))
  }
  out
}

s2c_long <- function(d, targets) {
  s2c_bind(lapply(s2c_metrics, function(metric) {
    x <- d; x$metric <- metric; x$value <- x[[metric]]
    x$target <- s2c_target(metric, x$n_chains, targets)
    # A discrete/bounded parameter can have undefined tail ESS while R-hat,
    # bulk ESS and mean MCSE remain valid. Assess each metric separately.
    invalid_draws <- x$unavailable_reason %in%
      c("nonfinite draws", "too few chains/draws", "constant sampled draws")
    x$assessable <- !invalid_draws & is.finite(x$value)
    x$departure <- s2c_departure(metric, x$value, x$target)
    x$shortfall_or_excess <- pmax(0, x$departure)
    x$rhat_minus_one <- if (metric == "rhat") x$value - 1 else NA_real_
    x$failed <- ifelse(x$assessable, s2c_failure(metric, x$value, x$target), NA)
    x
  }))
}
s2c_summary <- function(long) {
  all <- long; all$group <- paste0("Family: ", all$family)
  beta <- long[long$family == "Beta", , drop = FALSE]
  each <- beta; each$group <- paste0("Beta: ", each$predictor)
  slopes <- beta[beta$predictor != "(Intercept)", , drop = FALSE]; slopes$group <- "Beta: all slopes"
  x <- s2c_bind(list(all, each, slopes))
  groups <- split(x, interaction(x$model, x$group, x$metric, drop = TRUE))
  s2c_bind(lapply(groups, function(z) {
    v <- z$value[z$assessable]; bad <- z$failed[z$assessable]
    adverse_prob <- if (z$metric[1] %in% c("ess_bulk", "ess_tail")) .05 else .95
    base <- z[1, c("model", "atlas", "label", "variant", "group", "metric"), drop = FALSE]
    cbind(base, data.frame(n_parameters = nrow(z), n_assessed = length(v),
      n_unavailable = sum(!z$assessable), n_failed = sum(bad),
      fraction_failed_assessed = if (length(v)) mean(bad) else NA_real_,
      fraction_failed_all = sum(bad) / nrow(z),
      median = if (length(v)) median(v) else NA_real_,
      adverse_quantile = if (length(v)) unname(quantile(v, adverse_prob)) else NA_real_,
      worst = if (length(v)) if (adverse_prob == .05) min(v) else max(v) else NA_real_,
      median_departure = if (length(v)) median(z$departure[z$assessable]) else NA_real_,
      mean_shortfall_or_excess = if (length(v)) mean(z$shortfall_or_excess[z$assessable]) else NA_real_,
      target = z$target[1]))
  }))
}
s2c_compare_metadata <- function(metadata) {
  rows <- list()
  for (atlas in unique(vapply(metadata, `[[`, "", "atlas"))) {
    ms <- Filter(function(m) m$atlas == atlas, metadata)
    if (length(ms) < 2) next
    for (pair in combn(seq_along(ms), 2, simplify = FALSE)) {
      a <- ms[[pair[1]]]; b <- ms[[pair[2]]]
      sites <- setequal(a$sites, b$sites); species <- setequal(a$species, b$species)
      predictors <- setequal(a$predictors, b$predictors)
      x_equal <- sites && predictors && isTRUE(all.equal(unname(a$X),
        unname(b$X[match(a$sites, b$sites), match(a$predictors, b$predictors), drop = FALSE]), check.attributes = FALSE))
      y_equal <- sites && species && isTRUE(all.equal(unname(a$Y),
        unname(b$Y[match(a$sites, b$sites), match(a$species, b$species), drop = FALSE]), check.attributes = FALSE))
      rows[[length(rows) + 1L]] <- data.frame(atlas = atlas, model_a = a$model, model_b = b$model,
        same_cells = sites, same_species = species, same_predictors = predictors,
        same_design_values = x_equal, same_response_values = y_equal,
        same_scaling = isTRUE(all.equal(a$XScalePar, b$XScalePar)),
        same_sampling = identical(a$sampling, b$sampling),
        shared_cells = length(intersect(a$sites, b$sites)),
        cells_a = length(a$sites), cells_b = length(b$sites),
        note = "All spatial diagnostics use each model's own fitted cells.")
    }
  }
  s2c_bind(rows)
}
s2c_matched <- function(d, metadata) {
  out <- list()
  for (atlas in unique(d$atlas)) {
    models <- vapply(Filter(function(m) m$atlas == atlas, metadata), `[[`, "", "model")
    if (length(models) < 2) next
    for (pair in combn(models, 2, simplify = FALSE)) {
      a <- d[d$model == pair[1] & d$family != "Alpha (saved basis)", , drop = FALSE]
      b <- d[d$model == pair[2] & d$family != "Alpha (saved basis)", , drop = FALSE]
      keys <- c("family", "parameter", "predictor", "species", "level")
      # Blank optional keys may be read back from CSV as logical NA.
      for (v in c("predictor", "species", "level")) {
        a[[v]] <- as.character(a[[v]]); b[[v]] <- as.character(b[[v]])
        a[[v]][is.na(a[[v]])] <- ""; b[[v]][is.na(b[[v]])] <- ""
      }
      keep <- c(keys, "rhat", "ess_bulk", "ess_tail", "mcse_sd_ratio", "unavailable_reason")
      z <- merge(a[, keep], b[, keep], by = keys, all = TRUE, suffixes = c("_a", "_b"))
      z$atlas <- atlas; z$model_a <- pair[1]; z$model_b <- pair[2]
      z$match_status <- ifelse(is.na(z$unavailable_reason_a), "only model B",
                               ifelse(is.na(z$unavailable_reason_b), "only model A", "matched"))
      for (metric in s2c_metrics) z[[paste0(metric, "_change_b_minus_a")]] <- z[[paste0(metric, "_b")]] - z[[paste0(metric, "_a")]]
      out[[length(out) + 1L]] <- z
    }
  }
  s2c_bind(out)
}

s2c_colour_score <- function(metric, value, target) {
  # -1 green, 0 white at target, +1 red. ESS has the opposite direction.
  n <- max(length(metric), length(value), length(target))
  metric <- rep_len(metric, n); value <- rep_len(value, n); target <- rep_len(target, n)
  score <- ifelse(metric %in% c("ess_bulk", "ess_tail"), (target - value) / target,
    ifelse(metric == "rhat", (pmax(value, 1) - target) / pmax(target - 1, .Machine$double.eps),
           (value - target) / target))
  pmax(-1, pmin(1, score))
}
s2c_comparison_heatmap <- function(summary, metadata, targets, atlas_labels, output_file) {
  variables <- unique(unlist(lapply(metadata, `[[`, "predictors")))
  names_nice <- c("(Intercept)" = "Intercept", tmean_breeding = "Temperature",
    prec_breeding = "Precipitation", hh = "Land-use heterogeneity",
    perc_urban = "Urban cover", perc_cropland = "Cropland cover",
    perc_pasture = "Pasture cover", perc_grass_shrub = "Grass/shrubland cover",
    perc_forest_deciduous = "Deciduous forest cover",
    perc_forest_evergreen = "Evergreen forest cover", perc_forest_mixed_other = "Mixed/other forest cover")
  variable_names <- ifelse(variables %in% names(names_nice), names_nice[variables], variables)
  groups <- c("Family: Beta", paste0("Beta: ", variables),
              paste0("Family: ", c("Gamma", "V", "rho", "Alpha (saved basis)", "Residual correlation")))
  groups <- groups[groups %in% summary$group]
  row_labels <- c("Family: Beta" = "All Beta coefficients",
    setNames(paste0("Beta: ", variable_names), paste0("Beta: ", variables)),
    "Family: Gamma" = "Gamma (traits)", "Family: V" = "V (coefficient covariance)",
    "Family: rho" = "Phylogenetic rho", "Family: Alpha (saved basis)" = "Spatial ranges (saved factors)",
    "Family: Residual correlation" = "Residual species correlations")
  models <- vapply(metadata, `[[`, "", "model")
  labels <- vapply(metadata, function(m) {
    atlas <- if (m$atlas %in% names(atlas_labels)) atlas_labels[[m$atlas]] else m$atlas
    variant <- if (m$variant == "Environment only") "Env. only" else
      if (startsWith(m$variant, "Space nf")) sub("Space nf", "Spatial nf", m$variant) else m$variant
    paste(variant, atlas, sep = "\n")
  }, "")
  duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicated_labels] <- paste(labels[duplicated_labels],
    sub("_Atlas.*", "", models[duplicated_labels]), sep = "\n")
  skeleton <- expand.grid(model = models, group = groups, metric = s2c_metrics, stringsAsFactors = FALSE)
  z <- merge(skeleton, summary, by = c("model", "group", "metric"), all.x = TRUE)
  z$status <- ifelse(is.na(z$n_parameters), "Inapplicable",
    ifelse(z$n_assessed == 0, "Unavailable", ifelse(z$n_failed > 0, "Some parameters flagged", "None flagged")))
  z$colour_score <- s2c_colour_score(z$metric, z$median, z$target)
  value <- ifelse(z$metric == "rhat", sprintf("%.3f", z$median),
    ifelse(z$metric == "mcse_sd_ratio", sprintf("%.1f%%", z$median * 100), sprintf("%.0f", z$median)))
  z$label <- ifelse(z$status == "Inapplicable", "—",
    ifelse(z$status == "Unavailable", "Unavailable", paste0(value, "\n",
      sprintf("%.1f%%", z$fraction_failed_assessed * 100),
      ifelse(z$n_unavailable > 0, paste0("; NA", z$n_unavailable), ""))))
  ess_target <- targets$ess_per_chain * metadata[[1]]$row$n_chains
  metric_labels <- c(rhat = sprintf("Rank-normalised R-hat  (< %s)", targets$rhat),
    ess_bulk = sprintf("Bulk ESS  (≥ %s)", ess_target),
    ess_tail = sprintf("Tail ESS  (≥ %s)", ess_target),
    mcse_sd_ratio = sprintf("MCSE of mean / posterior SD  (≤ %g%%)", targets$mcse_sd_ratio * 100))
  z$metric <- factor(z$metric, levels = s2c_metrics, labels = metric_labels)
  z$column <- factor(z$model, levels = models, labels = labels)
  z$row <- factor(z$group, levels = rev(groups), labels = unname(row_labels[rev(groups)]))
  p <- ggplot2::ggplot(z, ggplot2::aes(column, row, fill = colour_score)) +
    ggplot2::geom_tile(colour = "white", linewidth = .7) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8, colour = "#202020", lineheight = .95) +
    ggplot2::facet_wrap(~metric, ncol = 2) +
    ggplot2::scale_fill_gradient2(low = "#66AD83", mid = "white", high = "#D66C60",
      midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1),
      labels = c("Better", "Target", "Worse"), na.value = "#E5E5E5") +
    ggplot2::scale_x_discrete(position = "top", expand = ggplot2::expansion(add = .02)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = .02)) +
    ggplot2::labs(title = "MCMC convergence: model comparison",
      subtitle = "All selected models shown together. Each cell: median diagnostic across parameters, then percentage flagged.",
      caption = paste("Diagnostics are calculated separately for each coefficient across retained chains, then summarised; not diagnostics of community-mean coefficients or VP.",
        "Colour follows the displayed median: green better, white at target, red worse. Percentages failing are separate. Grey: unavailable or inapplicable.",
        "Spatial-range factor labels are conditional on the saved basis. Dates, input fingerprints and full denominators are recorded in comparison.rds and comparison-summary.csv.", sep = "\n"),
      x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
      axis.text.x = ggplot2::element_text(size = 9), axis.text.y = ggplot2::element_text(size = 10),
      panel.spacing = grid::unit(1.4, "lines"), strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9), plot.margin = ggplot2::margin(12, 15, 12, 12))
  variants <- vapply(metadata, `[[`, "", "variant")
  dividers <- which(head(variants, -1) != tail(variants, -1)) + .5
  if (length(dividers)) p <- p + ggplot2::geom_vline(xintercept = dividers, colour = "white", linewidth = 1.4)
  ggplot2::ggsave(output_file, p, width = max(14, length(models) * 1.8 + 3),
                 height = max(10.8, length(groups) * .65 + 3), dpi = 200, bg = "white", limitsize = FALSE)
  invisible(z)
}

run_multimodel_diagnostics <- function(patterns, root, samples, thin, chains, targets,
                                      output_root = file.path(root, "multimodel-diagnostics"),
                                      atlas_labels = c(Atlas1 = "1970s", Atlas2 = "1990s", Atlas3 = "2010s")) {
  for (pkg in c("Hmsc", "posterior", "ggplot2"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  if (!length(patterns) || anyNA(patterns) || any(!nzchar(patterns))) stop("Supply nonempty date patterns.")
  counts <- c(samples, thin, chains)
  if (any(!is.finite(counts)) || any(counts < 1) || any(counts != floor(counts))) stop("Counts must be positive integers.")
  if (!all(c("rhat", "ess_per_chain", "mcse_sd_ratio") %in% names(targets)) ||
      any(lengths(targets) != 1L) || any(!is.finite(unlist(targets))) ||
      targets$rhat < 1 || targets$ess_per_chain <= 0 || targets$mcse_sd_ratio <= 0)
    stop("Require rhat >= 1, ess_per_chain > 0 and mcse_sd_ratio > 0.")
  folders <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  folders <- sort(unique(folders[vapply(basename(folders), function(f)
    any(vapply(patterns, grepl, logical(1), x = f, fixed = TRUE)), logical(1))]))
  if (!length(folders)) stop("No model folders matched under ", root)
  filename <- sprintf("HPC_samples_%04d_thin_%02d_chains_%d.Rdata", samples, thin, chains)
  inputs <- file.path(folders, "Models", "Fitted", filename)
  if (any(!file.exists(inputs))) stop("Missing selected fitted inputs:\n", paste(inputs[!file.exists(inputs)], collapse = "\n"))
  hashes <- setNames(unname(tools::md5sum(inputs)), normalizePath(inputs))
  metadata <- lapply(inputs, function(input) {
    m <- s2c_read(input); s2c_validate(m, samples, thin, chains); s2c_metadata(m, input)
  })
  metadata <- metadata[order(vapply(metadata, `[[`, 0L, "factor_order"),
                             vapply(metadata, `[[`, "", "atlas"), vapply(metadata, `[[`, "", "model"))]
  # Reuse scalar diagnostics only when the complete selected input set has matching fingerprints.
  cache_files <- c(file.path(output_root, "comparison.rds"),
    list.files(output_root, pattern = "diagnostic-summaries\\.rds$", recursive = TRUE, full.names = TRUE))
  result <- NULL
  for (file in cache_files[file.exists(cache_files)]) {
    candidate <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!is.null(candidate$diagnostics) && setequal(names(candidate$manifest$inputs), names(hashes)) &&
        identical(unname(candidate$manifest$inputs[names(hashes)]), unname(hashes))) {
      result <- candidate; message("Reusing validated parameter diagnostics: ", file); break
    }
  }
  if (is.null(result)) {
    diagnostic_rows <- list(); spatial <- list(compensation = data.frame(), fields = data.frame(), dimensions = data.frame())
    for (meta in metadata) {
      message("Computing diagnostics: ", meta$model)
      m <- s2c_read(meta$input); a <- s2c_parameter_families(m, targets)
      sp <- s2c_spatial(m, a$beta, a$beta_key, targets)
      diagnostic_rows[[meta$model]] <- s2c_tag(s2c_bind(list(a$diagnostics, sp$diagnostics)), meta)
      for (v in names(spatial)) spatial[[v]] <- s2c_bind(list(spatial[[v]], s2c_tag(sp[[v]], meta)))
      rm(m, a, sp); gc(verbose = FALSE)
    }
    result <- list(diagnostics = s2c_bind(diagnostic_rows), spatial = spatial)
  }
  result$summary <- s2c_summary(s2c_long(result$diagnostics, targets))
  result$comparability <- s2c_compare_metadata(metadata)
  date_tag <- s2c_date_tag(vapply(metadata, `[[`, "", "model"))
  result$manifest <- list(created = as.character(Sys.time()), patterns = patterns, inputs = hashes, date_tag = date_tag,
    targets = targets, models = s2c_bind(lapply(metadata, `[[`, "row")),
    sampling = c(samples = samples, thin = thin, chains = chains),
    spatial_domain = "all fitted cells separately per model",
    versions = sapply(c("Hmsc", "posterior", "ggplot2"), function(p) as.character(utils::packageVersion(p))))
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(output_root)
  figure <- file.path(out, paste0(date_tag, "-multimodel-convergence.png"))
  result$figure_cells <- s2c_comparison_heatmap(result$summary, metadata, targets, atlas_labels, figure)
  write.csv(result$figure_cells, file.path(out, "comparison-summary.csv"), row.names = FALSE, na = "NA")
  write.csv(result$comparability, file.path(out, "model-comparability.csv"), row.names = FALSE)
  if (!identical(unname(tools::md5sum(names(hashes))), unname(hashes))) stop("Fitted inputs changed.")
  saveRDS(result, file.path(out, "comparison.rds"))
  message("Saved one comparison figure: ", figure)
  invisible(out)
}
