#!/usr/bin/env Rscript
# Combined convergence and in-sample fit entrypoint. Explicit folder pattern required.
.libPaths(c('~/Rlibs', .libPaths()))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_path <- if(length(script_arg)) sub('^--file=','',script_arg[1]) else 'scripts/2-hmsc/model-diagnostics/run-diagnostics.R'
helper_dir <- dirname(normalizePath(script_path))
source(file.path(helper_dir,'helpers.R'))
run_diagnostics(commandArgs(trailingOnly=TRUE),helper_dir)
