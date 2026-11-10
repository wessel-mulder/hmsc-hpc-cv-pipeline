#!/usr/bin/env Rscript
# Combined convergence and in-sample fit entrypoint. Explicit folder pattern required.
.libPaths(c('~/Rlibs', .libPaths()))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_path <- if(length(script_arg)) sub('^--file=','',script_arg[1]) else 'scripts/2-hmsc/S3_check_model_fit_TensorFlow.R'
helper_dir <- file.path(dirname(normalizePath(script_path)),'model-diagnostics')
source(file.path(helper_dir,'helpers.R'))
run_diagnostics(commandArgs(trailingOnly=TRUE),helper_dir)
