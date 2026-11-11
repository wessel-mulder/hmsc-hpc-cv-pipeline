# Run from project root:
# Rscript scripts/5-vis/vp-distributions/vp-wilcoxon-environment-only.R
# Exploratory paired signed-rank tests on exported species-level VP estimates.
# Not a posterior uncertainty analysis or an omnibus distribution-equality test.
model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-06_05-11-30"
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(tmean_breeding="Temperature",prec_breeding="Precipitation",
  hh="Land-use heterogeneity",perc_urban="Urban cover",perc_cropland="Cropland cover",
  perc_pasture="Pasture cover",perc_forest_deciduous="Deciduous forest cover",
  perc_forest_evergreen="Evergreen forest cover",perc_forest_mixed_other="Mixed/other forest cover",
  perc_grass_shrub="Grass/shrubland cover")
reference_species <- NULL; vp_tables <- list()
for(a in 1:3) {
  model <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage_EnvironmentOnly", model_stamp, a)
  raw <- read.csv(file.path(model_root,model,"Results",paste0(model,"parameter_estimates_VP_.csv")),
                   row.names=1,check.names=FALSE)
  stopifnot(setequal(rownames(raw),c(names(variable_labels),"TjurR2")),ncol(raw)==155L,
    !anyDuplicated(colnames(raw)),all(vapply(raw,is.numeric,logical(1))),all(is.finite(as.matrix(raw))))
  if(is.null(reference_species)) reference_species <- colnames(raw)
  stopifnot(setequal(colnames(raw),reference_species))
  raw <- raw[,reference_species,drop=FALSE]
  relative <- as.matrix(raw[names(variable_labels),,drop=FALSE])
  stopifnot(all(relative>=0),all(relative<=1),max(abs(colSums(relative)-1))<1e-8)
  tjur <- as.numeric(raw["TjurR2",])
  stopifnot(all(tjur>=0),all(tjur<=1))
  vp_tables[[a]] <- list(scaled=sweep(relative,2,tjur,"*")*100)
}
pairs <- list(c(1L,2L),c(2L,3L),c(1L,3L))
rows <- list(); paired_rows <- list(); row_i <- 0L
for(scale_name in "scaled") for(pair in pairs) for(variable in names(variable_labels)) {
  row_i <- row_i+1L
  earlier <- vp_tables[[pair[1]]][[scale_name]][variable,reference_species]
  later <- vp_tables[[pair[2]]][[scale_name]][variable,reference_species]
  stopifnot(identical(names(earlier),names(later)),all(is.finite(earlier)),all(is.finite(later)))
  delta <- later-earlier
  nonzero <- delta[delta!=0]
  abs_rank <- rank(signif(abs(nonzero),7),ties.method="average")
  statistic <- sum(abs_rank[nonzero>0])
  rank_biserial <- if(length(nonzero)) (2*statistic-sum(abs_rank))/sum(abs_rank) else NA_real_
  if(length(nonzero)) {
    test <- wilcox.test(later,earlier,paired=TRUE,alternative="two.sided",mu=0,
      exact=FALSE,correct=TRUE,conf.int=TRUE,conf.level=.95,tol.root=1e-8,digits.rank=7)
    stopifnot(abs(unname(test$statistic)-statistic)<1e-10)
    # Independently verify equivalent one-sample-difference and shuffled-ID tests.
    direct <- wilcox.test(delta,mu=0,alternative="two.sided",exact=FALSE,correct=TRUE,digits.rank=7)
    reversed <- later[rev(seq_along(later))]
    rematched <- reversed[match(names(earlier),names(reversed))]
    shuffled <- wilcox.test(rematched,earlier,paired=TRUE,exact=FALSE,correct=TRUE,digits.rank=7)
    stopifnot(abs(test$p.value-direct$p.value)<1e-12,abs(test$p.value-shuffled$p.value)<1e-12)
    p_raw <- test$p.value; hl <- unname(test$estimate)
    lower <- test$conf.int[1]; upper <- test$conf.int[2]
    note <- "normal approximation; continuity correction; tie adjustment"
  } else {
    p_raw <- hl <- lower <- upper <- NA_real_
    note <- "All paired differences zero; signed-rank inference undefined"
  }
  scale_sd <- sd(delta)
  skewness <- if(scale_sd>0) mean((delta-mean(delta))^3)/scale_sd^3 else NA_real_
  pair_label <- paste0(atlas_labels[pair[1]]," to ",atlas_labels[pair[2]])
  rows[[row_i]] <- data.frame(scale=scale_name,pair=pair_label,variable=variable,
    n_species=length(delta),n_nonzero=length(nonzero),n_positive=sum(delta>0),
    n_negative=sum(delta<0),n_zero=sum(delta==0),
    n_tied_absolute_differences=length(nonzero)-length(unique(signif(abs(nonzero),7))),
    mean_change_pp=mean(delta),median_change_pp=median(delta),
    hl_pseudomedian_change_pp=hl,hl_ci_lower_pp=lower,hl_ci_upper_pp=upper,
    rank_biserial=rank_biserial,V=statistic,p_raw=p_raw,
    difference_skewness=skewness,method_note=note)
  paired_rows[[row_i]] <- data.frame(scale=scale_name,pair=pair_label,variable=variable,
    species=reference_species,earlier_vp_percent=as.numeric(earlier),
    later_vp_percent=as.numeric(later),difference_pp=as.numeric(delta))
}
wilcoxon_results <- do.call(rbind,rows)
paired_species_values <- do.call(rbind,paired_rows)
stopifnot(nrow(wilcoxon_results)==30L,nrow(paired_species_values)==30L*155L,
          all(wilcoxon_results$n_positive+wilcoxon_results$n_negative+wilcoxon_results$n_zero==155L))
# Scaled-only family requested: ten predictors x three atlas pairs = 30 tests.
# Holm controls familywise error under arbitrary test dependence if the marginal
# p-values are valid. It cannot repair dependence among species within a test.
wilcoxon_results$p_holm <- p.adjust(wilcoxon_results$p_raw,method="holm",n=30L)
stopifnot(all(wilcoxon_results$p_holm>=wilcoxon_results$p_raw,na.rm=TRUE))
wilcoxon_summary <- do.call(rbind,lapply("scaled",function(scale_name) {
  do.call(rbind,lapply(vapply(pairs,function(p)paste0(atlas_labels[p[1]]," to ",atlas_labels[p[2]]),character(1)),
    function(pair_label) {
      d <- wilcoxon_results[wilcoxon_results$scale==scale_name & wilcoxon_results$pair==pair_label,]
      data.frame(scale=scale_name,pair=pair_label,n_tests=nrow(d),
        raw_p_below_005=sum(d$p_raw<.05,na.rm=TRUE),holm_p_below_005=sum(d$p_holm<.05,na.rm=TRUE))
  }))
}))
message("Completed 30 scaled-only exploratory paired signed-rank tests; Holm family = all 30 comparisons.")
print(wilcoxon_summary,row.names=FALSE)
print(wilcoxon_results[,c("scale","pair","variable","median_change_pp","p_raw","p_holm")],row.names=FALSE)
