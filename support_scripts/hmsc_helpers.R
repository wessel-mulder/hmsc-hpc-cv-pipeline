fix_hpc_alpha_samples <- function(temp, samples, verbose = FALSE) {
  if (!is.matrix(temp[[1]][[1]]$Alpha)) {
    if (verbose) message("Alpha is not a matrix; no fix required")
    return(temp)
  }

  if (verbose) message("Alpha is a matrix; converting rows to latent-factor vectors")
  for (i in seq_len(samples)) {
    temp_alpha_mat <- temp[[1]][[i]]$Alpha
    temp[[1]][[i]]$Alpha <- lapply(
      seq_len(nrow(temp_alpha_mat)),
      function(p) temp_alpha_mat[p, ]
    )
  }
  temp
}

ecospat_boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100,
                          PEplot = TRUE, rm.duplicate = TRUE,
                          method = "spearman") {
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    round(pi / ei, 10)
  }

  mini <- min(fit, obs)
  maxi <- max(fit, obs)

  if (length(nclass) == 1) {
    if (nclass == 0) {
      if (window.w == "default") window.w <- (max(fit) - min(fit)) / 10
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w) / res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else {
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini) / nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else {
    vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }

  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")
  f <- f[to.keep]

  if (length(f) < 2) {
    b <- NA
  } else {
    r <- seq_along(f)
    if (rm.duplicate) {
      r <- seq_along(f)[f != c(f[-1], TRUE)]
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)
  }

  HS <- apply(interval, 1, sum) / 2
  if (length(nclass) == 1 && nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1
  }
  HS <- HS[to.keep]

  if (PEplot) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio",
         col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }

  list(F.ratio = f, cor = round(b, 3), HS = HS)
}
