source("scripts/2-hmsc/multimodel-diagnostics/helpers.R")
source("scripts/2-hmsc/multimodel-diagnostics/spatial-factor-correlations.R")
set.seed(369)
targets <- list(rhat = 1.01, ess_per_chain = 100, mcse_sd_ratio = .05)
stopifnot(s2c_date_tag(c("2026-11-06_Atlas1", "2026-11-02_Atlas1")) == "2026-11-02-06",
          s2c_date_tag(c("2026-11-06_a", "2026-11-06_b")) == "2026-11-06",
          s2c_date_tag(c("2026-10-30_a", "2026-11-02_a")) == "2026-10-30_2026-11-02")
make_factor_model <- function(factors) {
  n <- 40L; nc <- 4L
  x <- cbind("(Intercept)" = 1, tmean_breeding = rnorm(30), perc_forest_deciduous = rnorm(30))
  spatial <- max(factors) > 0
  pi <- sample(1:30)
  m <- list(X = x, covNames = colnames(x), spNames = c("a", "b"), samples = n, thin = 1L,
    Y = matrix(0, 30, 2, dimnames = list(NULL, c("a", "b"))),
    studyDesign = data.frame(site = paste0("cell", 1:30)),
    ranLevels = if (spatial) list(site = list(sDim = 2)) else list(),
    Pi = if (spatial) matrix(pi, 30, 1) else matrix(integer(), 30, 0))
  m$postList <- lapply(seq_len(nc), function(k) lapply(seq_len(n), function(i) {
    nf <- rep_len(factors, n)[i]
    if (!spatial) return(list(Beta = matrix(rnorm(6), 3, 2)))
    eta <- cbind(x[, 2] + rnorm(30, sd = .08), x[, 3] + rnorm(30, sd = .08),
                 -x[, 2] + rnorm(30, sd = .08))[, seq_len(nf), drop = FALSE]
    raw <- matrix(NA_real_, 30, nf); raw[pi, ] <- eta
    list(Beta = matrix(rnorm(6), 3, 2), Eta = list(raw), Lambda = list(matrix(rnorm(nf * 2), nf, 2)))
  }))
  m
}
specs <- list(0L, 1L, 3L, c(1L, 3L))
rows <- list(); model_rows <- list()
for (i in seq_along(specs)) {
  m <- make_factor_model(specs[[i]])
  s2c_validate(m, 40L, 1L, 4L)
  id <- paste0("2026-11-12_Atlas1_test", i)
  result <- s2c_factor_summary(m, list(model = id, atlas = "Atlas1"), .95, targets)
  s <- result$summary
  rows[[i]] <- s; model_rows[[i]] <- data.frame(model = id, atlas = "Atlas1", factor_order = max(specs[[i]]))
  if (i == 1) stopifnot(nrow(s) == 3, all(s$availability == "no spatial factors"), all(is.na(s$saved_factor)))
  else {
    stopifnot(nrow(s) == max(specs[[i]]) * 3, all(is.na(s$median[s$predictor == "(Intercept)"])))
    warm <- s[s$predictor == "tmean_breeding" & s$saved_factor == 1, ]
    stopifnot(warm$median > .98, warm$lower > .95, warm$interval_excludes_zero)
    direct <- vapply(m$postList, function(chain) vapply(chain, function(draw)
      cor(m$X[, 2], draw$Eta[[1]][m$Pi[, 1], 1]), numeric(1)), numeric(40))
    stopifnot(abs(warm$median - median(direct)) < 1e-12)
    if (i == 4) {
      cold <- s[s$predictor == "tmean_breeding" & s$saved_factor == 3, ]
      stopifnot(cold$n_present == 80, cold$n_total == 160, cold$n_finite == 80,
                cold$median < -.98, cold$diagnostic_flag,
                cold$availability == "factor absent in some draws")
    }
  }
}
file <- tempfile("s2c-factor-layout-", fileext = ".png")
cells <- s2c_factor_heatmap(s2c_bind(rows), s2c_bind(model_rows), .95, c(Atlas1 = "1970s"), file)
stopifnot(file.exists(file), length(levels(cells$block)) == 4, length(unique(cells$column_id)) == 8)
no_factors <- s2c_factor_heatmap(rows[[1]], model_rows[[1]], .95, c(Atlas1 = "1970s"), file)
stopifnot(all(is.na(no_factors$median)), length(unique(no_factors$column_id)) == 1)
unlink(file)
cat("Zero/one/three/changing factors, posterior intervals, Pi mapping, date labels and flexible plot layout passed.\n")
