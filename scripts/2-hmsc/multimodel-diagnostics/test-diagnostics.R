# Run from repository root. Synthetic tests write only to a temporary directory.
source("scripts/2-hmsc/multimodel-diagnostics/helpers.R")
set.seed(369)
targets <- list(rhat = 1.01, ess_per_chain = 100, mcse_sd_ratio = .05)
z <- matrix(rnorm(8000), 2000, 4)
d <- s2c_diagnose(z, targets)
stopifnot(d$unavailable_reason == "", identical(d$rhat, as.numeric(posterior::rhat(z))),
          identical(d$ess_bulk, as.numeric(posterior::ess_bulk(z))),
          identical(d$ess_tail, as.numeric(posterior::ess_tail(z))))
shift <- z; shift[, 1] <- shift[, 1] + 3
stopifnot(s2c_diagnose(shift, targets)$rhat > 1.01)
ar <- replicate(4, as.numeric(arima.sim(list(ar = .99), n = 2000)))
stopifnot(s2c_diagnose(ar, targets)$ess_bulk < 400)
stopifnot(s2c_diagnose(matrix(1, 40, 4), targets)$unavailable_reason == "constant sampled draws")
invalid <- z; invalid[1, 1] <- NA
stopifnot(s2c_diagnose(invalid, targets)$unavailable_reason == "nonfinite draws")
stopifnot(s2c_diagnose(z[, 1, drop = FALSE], targets)$unavailable_reason == "too few chains/draws")
bounded <- matrix(rep(c(.90, .95, 1, 1, 1), 200), 250, 4)
bd <- s2c_diagnose(bounded, targets)
stopifnot(is.finite(bd$rhat), is.finite(bd$ess_bulk), is.na(bd$ess_tail))
bl <- s2c_long(cbind(s2c_key("rho", "rho[1]"), bd), targets)
stopifnot(identical(bl$assessable, c(TRUE, TRUE, FALSE, TRUE)))
stopifnot(s2c_failure("rhat", 1.01, 1.01), !s2c_failure("ess_bulk", 400, 400),
          !s2c_failure("mcse_sd_ratio", .05, .05),
          s2c_departure("ess_tail", 300, 400) == 100,
          s2c_target("ess_bulk", c(2, 4), targets)[2] == 400)
key <- s2c_key("Beta", "temperature | a", "temperature", "a")
long <- s2c_long(cbind(data.frame(model = "a", atlas = "Atlas1", label = "a", variant = "a"), key, d), targets)
strict <- targets; strict$rhat <- d$rhat - .00001
stopifnot(s2c_long(cbind(key, d), strict)$failed[1], !long$failed[1])
x <- cbind("(Intercept)" = 1, temperature = rnorm(20), forest = rnorm(20))
eta <- matrix(rnorm(40), 20, 2); lambda <- matrix(rnorm(6), 2, 3)
q <- qr.Q(qr(matrix(rnorm(4), 2)))
stopifnot(max(abs(s2c_spatial_cor(x[, 2:3], eta, lambda) - cor(x[, 2:3], eta %*% lambda))) < 1e-12,
          max(abs(s2c_spatial_cor(x[, 2:3], eta, lambda) - s2c_spatial_cor(x[, 2:3], eta %*% q, t(q) %*% lambda))) < 1e-12,
          max(abs(s2c_project(x, eta, lambda) - qr.coef(qr(x), eta %*% lambda))) < 1e-12,
          is.null(s2c_project(cbind(x, x[, 2]), eta, lambda)))
make_model <- function(spatial = FALSE, changing = FALSE) {
  m <- list(covNames = colnames(x), spNames = c("a", "b", "c"), X = x, XScalePar = rbind(0, 1),
    Y = matrix(rep(c(0, 1), 30), 20, 3, dimnames = list(NULL, c("a", "b", "c"))),
    studyDesign = data.frame(site = paste0("cell", 1:20)), samples = 40L, thin = 1L,
    distr = matrix(0, 3, 2, dimnames = list(NULL, c("family", "variance"))),
    ranLevels = if (spatial) list(site = list(sDim = 2, alphapw = cbind(c(0, 100), c(.5, .5)))) else list(),
    Pi = if (spatial) matrix(1:20, 20, 1) else matrix(integer(), 20, 0),
    rhopw = cbind(1, 1), trNames = c("(Intercept)", "trait"))
  m$postList <- lapply(1:4, function(k) lapply(1:40, function(i) {
    nf <- if (changing && i %% 2 == 0) 1L else 2L
    list(Beta = matrix(rnorm(9), 3, 3), Gamma = matrix(rnorm(6), 3, 2), rho = 1,
      sigma = rep(1, 3), Eta = if (spatial) list(matrix(rnorm(20 * nf), 20, nf)) else NULL,
      Lambda = if (spatial) list(matrix(rnorm(3 * nf), nf, 3)) else NULL,
      Alpha = if (spatial) list(sample(1:2, nf, replace = TRUE)) else NULL)
  }))
  m
}
m <- make_model(); a <- s2c_parameter_families(m, targets)
stopifnot(identical(as.numeric(a$beta[3, 2, ]), as.vector(m$postList[[2]][[3]]$Beta)),
          !any(a$diagnostics$family %in% c("rho", "sigma")),
          grepl("inapplicable", s2c_spatial(m, a$beta, a$beta_key, targets)$notes))
stopifnot(inherits(try(s2c_validate(m, 41L, 1L, 4L), silent = TRUE), "try-error"))
# Matrix subsetting drops model.matrix's assign attribute; values remain comparable.
metadata_input <- file.path(tempfile("s2c-metadata-"), "2026-11-12_Atlas1_A", "Models", "Fitted", "test.Rdata")
dir.create(dirname(metadata_input), recursive = TRUE)
writeLines("metadata fixture", metadata_input)
ma <- s2c_metadata(m, metadata_input)
mb <- ma; mb$model <- "other"; attr(ma$X, "assign") <- 0:2
mb$sites <- rev(ma$sites); mb$X <- ma$X[20:1, , drop = FALSE]; mb$Y <- ma$Y[20:1, , drop = FALSE]
stopifnot(s2c_compare_metadata(list(ma, mb))$same_design_values,
          s2c_compare_metadata(list(ma, mb))$same_response_values)
mb$X[1, 2] <- mb$X[1, 2] + 1
stopifnot(!s2c_compare_metadata(list(ma, mb))$same_design_values)
# Undefined diagnostics retain their denominator and never count as passes.
invalid_row <- cbind(data.frame(model = "a", atlas = "Atlas1", label = "a", variant = "a"),
                     key, s2c_diagnose(matrix(1, 40, 4), targets))
sm <- s2c_summary(s2c_long(invalid_row, targets))
stopifnot(all(sm$n_assessed == 0), all(sm$n_unavailable == 1), all(is.na(sm$fraction_failed_assessed)))
for (changing in c(FALSE, TRUE)) {
  m <- make_model(TRUE, changing); a <- s2c_parameter_families(m, targets)
  sp <- s2c_spatial(m, a$beta, a$beta_key, targets)
  stopifnot(nrow(sp$compensation) == 36, sum(sp$diagnostics$family == "Residual correlation") == 3,
            sum(sp$diagnostics$family == "Alpha (saved basis)") == if (changing) 0 else 2)
  if (changing) stopifnot(length(unique(sp$dimensions$n_factors)) == 2)
  else stopifnot(all(sp$diagnostics$median[sp$diagnostics$family == "Alpha (saved basis)"] %in% c(0, 50, 100)))
  m$X[, 3] <- m$X[, 2]
  sp <- s2c_spatial(m, a$beta, a$beta_key, targets)
  stopifnot(nrow(sp$compensation) == 0, any(grepl("Rank-deficient", sp$notes)))
}
m <- make_model(TRUE); m$ranLevels[[1]]$alphapw[, 2] <- c(1, 0)
for (k in seq_along(m$postList)) for (i in seq_along(m$postList[[k]])) m$postList[[k]][[i]]$Alpha[[1]][] <- 1L
a <- s2c_parameter_families(m, targets); sp <- s2c_spatial(m, a$beta, a$beta_key, targets)
stopifnot(!any(sp$diagnostics$family == "Alpha (saved basis)"))
cat("Direct diagnostic, invalid-draw, benchmark, extraction, spatial invariance, changing-factor and rank-deficiency checks passed.\n")
if ("--integration" %in% commandArgs(trailingOnly = TRUE)) {
  root <- tempfile("s2c-fixture-"); dir.create(root)
  for (spatial in c(FALSE, TRUE)) {
    folder <- file.path(root, paste0("2026-11-12_Atlas1_", if (spatial) "Space" else "EnvironmentOnly"), "Models", "Fitted")
    dir.create(folder, recursive = TRUE)
    fitted_model <- list(posteriors = make_model(spatial))
    save(fitted_model, file = file.path(folder, "HPC_samples_0040_thin_01_chains_4.Rdata"))
  }
  out <- run_multimodel_diagnostics(c("2026-11-12", "2026-11-12"), root, 40L, 1L, 4L, targets)
  result <- readRDS(file.path(out, "comparison.rds"))
  stopifnot(nrow(result$manifest$models) == 2, nrow(result$comparability) == 1,
            !any(result$diagnostics$family == "sigma"),
            all(result$manifest$inputs == unname(tools::md5sum(names(result$manifest$inputs)))),
            length(list.files(out, pattern = "\\.png$", recursive = TRUE)) == 1,
            length(unique(result$figure_cells$column)) == 2,
            length(unique(result$figure_cells$metric)) == 4)
  meta <- lapply(seq_len(nrow(result$manifest$models)), function(i) as.list(result$manifest$models[i, ]))
  matched <- s2c_matched(result$diagnostics, meta)
  stopifnot(sum(matched$family == "Beta" & matched$match_status == "matched") == 9)
  relaxed <- targets; relaxed$rhat <- 1.1
  run_multimodel_diagnostics("2026-11-12", root, 40L, 1L, 4L, relaxed)
  updated <- readRDS(file.path(out, "comparison.rds"))
  stopifnot(identical(result$diagnostics, updated$diagnostics), updated$manifest$targets$rhat == 1.1,
            length(list.files(out, pattern = "\\.png$", recursive = TRUE)) == 1,
            sum(updated$summary$n_failed[updated$summary$metric == "rhat"]) <=
              sum(result$summary$n_failed[result$summary$metric == "rhat"]))
  cat("Full two-model integration passed; output:", out, "\n")
}
args <- commandArgs(trailingOnly = TRUE)
if ("--models" %in% args) {
  i <- match("--models", args)
  if (i == length(args)) stop("Supply the completed output directory after --models.")
  output <- args[i + 1L]; result <- readRDS(file.path(output, "comparison.rds"))
  manifest <- result$manifest; d <- result$diagnostics
  stopifnot(identical(unname(tools::md5sum(names(manifest$inputs))), unname(manifest$inputs)))
  for (j in seq_len(nrow(manifest$models))) {
    meta <- manifest$models[j, ]; m <- s2c_read(meta$input)
    z <- d[d$model == meta$model, ]
    # Compare actual extraction and diagnostics for every output family.
    for (family in unique(z$family)) {
      selected_row <- z[which(z$family == family)[1], ]
      p <- match(selected_row$predictor, m$covNames); s <- match(selected_row$species, m$spNames)
      level <- match(selected_row$level, names(m$ranLevels))
      getter <- switch(family,
        Beta = function(draw) draw$Beta[p, s],
        Gamma = function(draw) draw$Gamma[1, 1],
        V = function(draw) draw$V[1, 1],
        rho = function(draw) draw$rho[1],
        sigma = function(draw) draw$sigma[match(selected_row$species, m$spNames)],
        `Alpha (saved basis)` = function(draw) m$ranLevels[[level]]$alphapw[draw$Alpha[[level]][1], 1],
        `Residual correlation` = function(draw) {
          covariance <- crossprod(draw$Lambda[[level]])
          covariance[1, 2] / sqrt(covariance[1, 1] * covariance[2, 2])
        },
        `Spatial overlap` = function(draw) {
          effect <- draw$Eta[[level]][m$Pi[, level], , drop = FALSE] %*% draw$Lambda[[level]]
          cor(m$X[, p], effect[, s])
        },
        `Projected space` = function(draw) {
          effect <- draw$Eta[[level]][m$Pi[, level], , drop = FALSE] %*% draw$Lambda[[level]]
          qr.coef(qr(m$X), effect)[p, s]
        },
        `Beta + projected space` = function(draw) {
          effect <- draw$Eta[[level]][m$Pi[, level], , drop = FALSE] %*% draw$Lambda[[level]]
          draw$Beta[p, s] + qr.coef(qr(m$X), effect)[p, s]
        })
      stopifnot(!is.null(getter))
      draws <- vapply(m$postList, function(chain) vapply(chain, getter, numeric(1)), numeric(m$samples))
      if (family %in% c("Projected space", "Beta + projected space", "Spatial overlap")) {
        # First establish agreement with the independent dense-effect calculation.
        # Then retain exact low-rank arithmetic for rank diagnostics: floating-point
        # differences at machine precision can break otherwise tied ranks.
        lowrank_getter <- function(draw) {
          eta <- draw$Eta[[level]][m$Pi[, level], , drop = FALSE]; lambda <- draw$Lambda[[level]]
          if (family == "Spatial overlap") {
            variable <- which(apply(m$X, 2, sd) > 0)
            return(s2c_spatial_cor(m$X[, variable, drop = FALSE], eta, lambda)[match(p, variable), s])
          }
          projected <- (qr.coef(qr(m$X), eta) %*% lambda)[p, s]
          if (family == "Projected space") projected else draw$Beta[p, s] + projected
        }
        exact <- vapply(m$postList, function(chain) vapply(chain, lowrank_getter, numeric(1)), numeric(m$samples))
        stopifnot(max(abs(draws - exact), na.rm = TRUE) < 1e-12)
        draws <- exact
      }
      direct <- s2c_diagnose(draws, manifest$targets)
      stopifnot(isTRUE(all.equal(direct$rhat, as.numeric(posterior::rhat(draws)))),
                isTRUE(all.equal(direct$ess_bulk, as.numeric(posterior::ess_bulk(draws)))))
      # Projection/correlation low-rank algebra can change ranks of nearly equal values very slightly.
      for (v in c("median", "sd", "lower", "upper", s2c_metrics, "mcse_lower", "mcse_upper"))
        if (!isTRUE(all.equal(selected_row[[v]], direct[[v]], tolerance = 1e-7)))
          stop(meta$model, ": ", family, " ", v, " saved=", selected_row[[v]], " direct=", direct[[v]])
    }
    message("Direct saved-draw verification passed: ", meta$model)
    rm(m); gc(verbose = FALSE)
  }
  metadata <- lapply(seq_len(nrow(manifest$models)), function(i) as.list(manifest$models[i, ]))
  matched <- s2c_matched(d, metadata)
  pairs <- unique(matched[, c("model_a", "model_b")])
  for (i in seq_len(nrow(pairs))) {
    a <- d$parameter[d$model == pairs$model_a[i] & d$family == "Beta"]
    b <- d$parameter[d$model == pairs$model_b[i] & d$family == "Beta"]
    beta <- matched[matched$model_a == pairs$model_a[i] & matched$model_b == pairs$model_b[i] & matched$family == "Beta", ]
    stopifnot(sum(beta$match_status == "matched") == length(intersect(a, b)),
              nrow(beta) == length(union(a, b)))
  }
  stopifnot(identical(unname(tools::md5sum(names(manifest$inputs))), unname(manifest$inputs)))
  cat("Real-fit extraction, direct family diagnostics, matching and unchanged fingerprints passed.\n")
}
