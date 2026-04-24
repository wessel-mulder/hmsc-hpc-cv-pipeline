model_fit_limits <- function(type, cMF, cMFCV) {
  switch(
    type,
    "RMSE" = c(rep(max(cMF[[type]], cMFCV[[type]], na.rm = TRUE), 2), 0, 0),
    "AUC" = c(0, 1, 0.5, 0.5),
    "O.AUC" = c(0, 1, 0.5, 0.5),
    c(rep(1, 2), 0, 0)
  )
}

plot_model_fit_cv <- function(type, cMF, cMFCV, modelnames, thin, samples) {
  limit <- model_fit_limits(type, cMF, cMFCV)
  plot(
    cMF[[type]], cMFCV[[type]],
    xlim = c(-limit[1], limit[2]),
    ylim = c(-limit[1], limit[2]),
    xlab = "explanatory power",
    ylab = "predictive power",
    main = sprintf(
      "%s:\n%s thin = %i, samples = %i\nMF: mean(%1$s) = %.4f MFSCV: mean(%1$s) = %.4f",
      type, modelnames, thin, samples,
      mean(cMF[[type]], na.rm = TRUE),
      mean(cMFCV[[type]], na.rm = TRUE)
    )
  )
  abline(0, 1)
  abline(v = limit[3])
  abline(h = limit[4])
}
