# Explore land coverage among atlas cells with good survey coverage.
# Run from the repository root: Rscript assessing-ocean-threshold/assess-ocean-threshold.R
suppressPackageStartupMessages({library(readxl); library(ggplot2); library(dplyr)})

root <- 'data-new/0-pre'
env_path <- file.path(root, 'derived/environmental/X_environment.csv')
effort_path <- file.path(root, 'raw/effort/dækning atlas1+2.xlsx')
out_dir <- 'assessing-ocean-threshold'
fig_dir <- file.path(out_dir, 'figures')
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

env <- read.csv(env_path, check.names=FALSE)
required <- c('survey','kvadratkod','atlas','pct_class_0','pct_class_77')
stopifnot(all(required %in% names(env)), !anyDuplicated(env$survey))

# Good effort only: Atlas 1 uses UNDS_KODE_; Atlas 2 uses corrected UNDS_KOD_1.
a1 <- read_excel(effort_path, sheet='kvadrater atlas1')
a2 <- read_excel(effort_path, sheet='kvadrater atlas2')
good_surveys <- c(paste0(a1[['KVADRATKOD']][a1[['UNDS_KODE_']]==1], '_1'),
                  paste0(a2[['KVADRATKOD']][a2[['UNDS_KOD_1']]==1], '_2'))
env$coverage <- ifelse(env$atlas == 3, 'Assumed good (Atlas 3)',
                       ifelse(env$survey %in% good_surveys, 'Good', 'Not good'))
env$water_pct <- env$pct_class_0 + env$pct_class_77
env$land_pct <- 100 - env$water_pct

# Keep only measured Good Atlas 1/2 cells and all Atlas 3 cells under the
# same explicit assumed-good convention used by preprocessing.
d <- env[env$coverage != 'Not good' & is.finite(env$land_pct) &
           env$land_pct >= 0 & env$land_pct <= 100, ]
d$atlas_label <- factor(paste0('Atlas ', d$atlas), levels=paste0('Atlas ', 1:3))

thresholds <- seq(0, 50, by=1)
threshold_summary <- do.call(rbind, lapply(split(d, d$atlas_label), function(z) {
  data.frame(atlas=as.character(z$atlas_label[1]), threshold_land_pct=thresholds,
             n_cells=nrow(z), n_retained=sapply(thresholds, function(t) sum(z$land_pct >= t)),
             pct_retained=sapply(thresholds, function(t) 100*mean(z$land_pct >= t)))
}))
write.csv(threshold_summary, file.path(out_dir, 'threshold-retention-summary.csv'), row.names=FALSE)
write.csv(d[,c('survey','kvadratkod','atlas','coverage','pct_class_0','pct_class_77','water_pct','land_pct')],
          file.path(out_dir, 'good-coverage-cell-land-cover.csv'), row.names=FALSE)

stats <- do.call(rbind, lapply(split(d, d$atlas_label), function(z) {
  data.frame(atlas=as.character(z$atlas_label[1]), coverage=z$coverage[1], n_cells=nrow(z),
             median_land_pct=median(z$land_pct), mean_land_pct=mean(z$land_pct),
             below_25_land=sum(z$land_pct < 25), exactly_25_land=sum(z$land_pct == 25),
             at_least_25_land=sum(z$land_pct >= 25),
             pct_at_least_25_land=100*mean(z$land_pct >= 25))
}))
write.csv(stats, file.path(out_dir, 'land-coverage-summary.csv'), row.names=FALSE)

p1 <- ggplot(d, aes(x=land_pct, fill=atlas_label, colour=atlas_label)) +
  geom_density(alpha=.18, linewidth=.8, adjust=1.1) +
  geom_vline(xintercept=25, linetype='dashed', colour='#b23a48', linewidth=.7) +
  annotate('text', x=26, y=Inf, label='Current minimum: 25% land', hjust=0, vjust=1.3,
           colour='#9f2938', size=3.4) +
  scale_x_continuous(limits=c(0,100), breaks=seq(0,100,20), expand=c(0,0)) +
  scale_fill_brewer(palette='Dark2') + scale_colour_brewer(palette='Dark2') +
  labs(title='Land coverage in well-covered atlas cells',
       subtitle='Land = 100% − ocean (class 0) − freshwater (class 77)',
       x='Land coverage per cell (%)', y='Density', colour='Atlas', fill='Atlas',
       caption='Atlas 1–2: Good survey coverage only. Atlas 3: coverage assumed good, as in preprocessing. Dashed line marks the inclusive 25% land cutoff.') +
  theme_minimal(base_size=12) + theme(legend.position='top', panel.grid.minor=element_blank())
ggsave(file.path(fig_dir, 'land-coverage-distribution.png'), p1, width=9, height=5.5, dpi=300)

p2 <- ggplot(threshold_summary, aes(threshold_land_pct, pct_retained, colour=atlas)) +
  geom_line(linewidth=1) + geom_point(data=subset(threshold_summary, threshold_land_pct %in% c(10,20,25,30,40)), size=1.5) +
  geom_vline(xintercept=25, linetype='dashed', colour='#555555') +
  scale_x_continuous(breaks=seq(0,50,5)) + scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
  scale_colour_brewer(palette='Dark2') +
  labs(title='Cells retained under alternative minimum-land thresholds',
       x='Minimum land coverage (%)', y='Cells retained (%)', colour='Atlas',
       caption='Retention uses land coverage ≥ threshold; at 25%, this matches the current rule of ocean + freshwater ≤ 75%.') +
  theme_minimal(base_size=12) + theme(legend.position='top', panel.grid.minor=element_blank())
ggsave(file.path(fig_dir, 'threshold-sensitivity.png'), p2, width=9, height=5.5, dpi=300)

cat('Good-coverage cells by atlas and 25% land retention:\n')
print(stats, row.names=FALSE)
