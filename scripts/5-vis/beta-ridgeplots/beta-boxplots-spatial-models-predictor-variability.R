# Run from project root:
# Rscript scripts/5-vis/beta-ridgeplots/beta-boxplots-spatial-models-predictor-variability.R
# Betas retain the fitted atlas-specific predictor scaling: effect per one
# within-atlas predictor SD, on the probit linear-predictor scale.
# No conversion to original predictor units and no significance filtering.

required <- c("ggplot2", "ggridges", "posterior")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(ggridges)})
output_dir <- "scripts/5-vis/beta-ridgeplots"
model_root <- "data-new/2-hmsc/derived/models"
model_stamp <- "2026-11-02_06-55-21"
atlas_labels <- c("1970s", "1990s", "2010s")
variable_labels <- c(
  tmean_breeding = "Temperature", prec_breeding = "Precipitation",
  hh = "Land-use heterogeneity", perc_urban = "Urban cover", perc_cropland = "Cropland cover",
  perc_pasture = "Pasture cover", perc_forest_deciduous = "Deciduous forest cover",
  perc_forest_evergreen = "Evergreen forest cover", perc_forest_mixed_other = "Mixed/other forest cover",
  perc_grass_shrub = "Grass/shrubland cover")
variable_colours <- c(
  "Temperature" = "firebrick3", "Precipitation" = "dodgerblue3",
  "Land-use heterogeneity" = "orchid3", "Urban cover" = "snow3",
  "Cropland cover" = "goldenrod1", "Pasture cover" = "darkorange",
  "Deciduous forest cover" = "#A6D96A", "Evergreen forest cover" = "#006400",
  "Mixed/other forest cover" = "#87CEEB", "Grass/shrubland cover" = "#B39DDB")
input_env <- new.env(parent = emptyenv())
load("data-new/2-hmsc/derived/preprocessed_data_min_occs_is_5_coverage_is_good_and_average_forest_split.RData",
     envir = input_env)
beta_draws <- list(); mean_rows <- list(); audit_rows <- list(); diagnostic_rows <- list()
reference_species <- NULL
for (atlas in 1:3) {
  folder <- sprintf("%s_Atlas%d_MinOccs5_CoverageGoodAverage", model_stamp, atlas)
  e <- new.env(parent = emptyenv())
  load(file.path(model_root, folder, "Models", "Fitted", "HPC_samples_0250_thin_100_chains_4.Rdata"), e)
  m <- e$fitted_model$posteriors
  predictors <- setdiff(m$covNames, "(Intercept)")
  stopifnot(m$nr == 1L, m$ns == 155L, setequal(predictors, names(variable_labels)),
            length(m$postList) == 4L, all(lengths(m$postList) == 250L))
  if (is.null(reference_species)) reference_species <- m$spNames
  stopifnot(setequal(m$spNames, reference_species), !anyDuplicated(m$spNames),
            all(rownames(m$XData) %in% rownames(input_env$X)),
            identical(rownames(m$XData), rownames(m$X)))
  pred_ix <- match(names(variable_labels), m$covNames)
  species_ix <- match(reference_species, m$spNames)
  raw <- input_env$X[rownames(m$XData), names(variable_labels), drop = FALSE]
  fitted_x <- as.matrix(m$X[, pred_ix, drop = FALSE])
  scaled_raw <- scale(raw)
  max_error <- max(abs(scaled_raw - fitted_x))
  stopifnot(is.finite(max_error), max_error < 1e-8,
            max(abs(colMeans(fitted_x))) < 1e-8,
            max(abs(apply(fitted_x, 2, sd) - 1)) < 1e-8)
  audit_rows[[atlas]] <- data.frame(atlas = atlas, variable = names(variable_labels),
    n_cells = nrow(raw), raw_mean = colMeans(raw), raw_sd = vapply(raw, sd, numeric(1)),
    fitted_mean = colMeans(fitted_x), fitted_sd = apply(fitted_x, 2, sd), max_scaling_error = max_error)
  draws <- array(NA_real_, c(250L, 4L, 10L, 155L),
    dimnames = list(NULL, NULL, names(variable_labels), reference_species))
  for (chain in 1:4) for (iteration in 1:250) {
    b <- m$postList[[chain]][[iteration]]$Beta
    stopifnot(identical(dim(b), c(11L, 155L)), all(is.finite(b)))
    draws[iteration, chain, , ] <- b[pred_ix, species_ix, drop = FALSE]
  }
  beta_draws[[atlas]] <- draws
  means <- apply(draws, c(3, 4), mean)
  # Verify selected species/predictor means by direct original-chain extraction.
  for (v in c(1L, 2L, 10L)) for (s in c(1L, 78L, 155L)) {
    direct <- unlist(lapply(m$postList, function(ch) vapply(ch,
      function(d) d$Beta[pred_ix[v], species_ix[s]], numeric(1))))
    stopifnot(abs(mean(direct) - means[v, s]) < 1e-12)
  }
  mean_rows[[atlas]] <- data.frame(atlas = atlas,
    variable = rep(names(variable_labels), times = 155L),
    species = rep(reference_species, each = 10L), beta_mean = as.vector(means))
  diagnostic_rows[[atlas]] <- do.call(rbind, lapply(1:10, function(v) {
    rh <- vapply(1:155, function(s) posterior::rhat(draws[, , v, s]), numeric(1))
    stopifnot(all(is.finite(rh)))
    data.frame(atlas = atlas, variable = names(variable_labels)[v],
      median_rhat = median(rh), max_rhat = max(rh), n_rhat_ge_1_1 = sum(rh >= 1.1))
  }))
  message("Atlas ", atlas, ": 155 species, ten predictors, 1,000 draws; scaling error ", signif(max_error, 3))
}
species_beta_means <- do.call(rbind, mean_rows)
predictor_scale_audit <- do.call(rbind, audit_rows)
beta_diagnostic_summary <- do.call(rbind, diagnostic_rows)
stopifnot(nrow(species_beta_means) == 4650L,
          all(table(species_beta_means$atlas, species_beta_means$variable) == 155L))


# One species-level posterior mean per box observation. Retain all draws above
# for diagnostics; do not pool posterior draws as independent observations.
box_data <- species_beta_means
box_data$variable <- factor(box_data$variable,levels=names(variable_labels),
                            labels=unname(variable_labels))
box_data$atlas <- factor(box_data$atlas,levels=1:3,labels=atlas_labels)
box_data$x_position <- as.integer(box_data$variable)*4 +
  c(-1.3,0,1.3)[as.integer(box_data$atlas)]
box_data$box_group <- interaction(box_data$variable,box_data$atlas,drop=TRUE)
# SDs from raw predictors in each atlas's actual fitted training sample,
# matching the physical SD used for each model coefficient's unit.
variability <- predictor_scale_audit
baseline_sd <- setNames(variability$raw_sd[variability$atlas==1],
                        variability$variable[variability$atlas==1])
variability$sd_ratio <- variability$raw_sd/baseline_sd[variability$variable]
stopifnot(all(is.finite(variability$sd_ratio)),
          all(variability$sd_ratio[variability$atlas==1]==1))
variability$variable <- factor(variability$variable,levels=names(variable_labels),
                               labels=unname(variable_labels))
variability$x_position <- as.integer(variability$variable)*4 +
  c(-1.3,0,1.3)[variability$atlas]
# Display transform only: ratio 1 aligns with Beta 0. The right-hand axis
# reverses this transform; ratios are not coefficient estimates.
variability$plot_y <- variability$sd_ratio-1
p <- ggplot(box_data,aes(x=x_position,y=beta_mean,fill=variable,alpha=atlas,group=box_group)) +
  geom_hline(yintercept=0,linetype='dashed',colour='grey45',linewidth=.4) +
  geom_boxplot(orientation='x',width=.65,colour='grey35',linewidth=.4,
               outlier.size=.8,outlier.alpha=.35) +
  stat_summary(fun=mean,geom='point',orientation='x',shape=23,size=1.8,
                fill='white',colour='grey25',alpha=1) +
  geom_line(data=variability,aes(x=x_position,y=plot_y,group=variable),
    inherit.aes=FALSE,colour='black',linewidth=.65) +
  geom_point(data=variability,aes(x=x_position,y=plot_y),inherit.aes=FALSE,
    shape=21,fill='white',colour='black',size=2.3,stroke=.7) +
  scale_fill_manual(values=variable_colours,guide='none') +
  scale_alpha_manual(name='Atlas (left to right)',values=c('1970s'=.45,'1990s'=.7,'2010s'=1)) +
  scale_x_continuous(breaks=seq_along(variable_labels)*4,
    labels=vapply(unname(variable_labels),function(z)paste(strwrap(z,width=17),collapse='\n'),character(1)),
    expand=expansion(add=1.8)) +
  scale_y_continuous(expand=expansion(mult=c(.04,.06)),
    sec.axis=sec_axis(~ . + 1,name='Predictor SD / Atlas 1 predictor SD',
      breaks=seq(0,3,by=.25),labels=function(x)paste0(format(x,trim=TRUE),'×'))) +
  labs(title='Spatial models: species–environment effect sizes',
    subtitle='155 species per box · Black lines and open circles: relative predictor variability (right axis)',
    x=NULL,y='Beta per 1 within-atlas predictor SD (probit scale)',
    caption=paste0('November 2 spatial models; all four chains × 250 retained draws. Within each variable: 1970s, 1990s, 2010s.\n',
      'Boxes: median and IQR; whiskers: 1.5 × IQR; dots: outliers; white diamonds: means; dashed line: zero.\n',
      'Boxes describe variation among species’ posterior mean coefficients, not posterior uncertainty. No coefficient filtering.\n',
      'Black lines: raw predictor SD in each atlas’s fitted training cells / Atlas 1 SD; 1× means unchanged variability (right axis).\n',
      'SD ratios are descriptive and use each atlas’s training sample. Lines connect atlas values only within each predictor.')) +
  theme_minimal(base_size=12) + theme(legend.position='top',legend.justification='left',
    panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
    axis.text.x=element_text(size=10),plot.title=element_text(face='bold',size=17),
    plot.caption=element_text(hjust=0,size=9),plot.margin=margin(12,12,12,12))
boxes <- ggplot_build(p)$data[[2]]
expected <- as.numeric(tapply(box_data$beta_mean,list(box_data$variable,box_data$atlas),median))
stopifnot(nrow(boxes)==30L,max(abs(sort(boxes$middle)-sort(expected)))<1e-10)
path <- file.path(output_dir,'beta-boxplots-spatial-models-predictor-variability.png')
ggsave(path,p,width=18,height=8,dpi=300,bg='white')
print(variability[,c('atlas','variable','n_cells','raw_sd','sd_ratio')],row.names=FALSE)
message('Saved ',path,'; all 30 medians checked; no original ridgeplots overwritten.')
