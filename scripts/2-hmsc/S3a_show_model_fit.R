#!/usr/bin/env Rscript
# Compatibility entrypoint: identical positional arguments and MF_preCV cache format.
.libPaths(c('~/Rlibs', .libPaths()))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
script_path <- if(length(script_arg)) sub('^--file=','',script_arg[1]) else 'scripts/2-hmsc/S3a_show_model_fit.R'
helper_dir <- file.path(dirname(normalizePath(script_path)),'model-diagnostics')
source(file.path(helper_dir,'helpers.R'))
run_diagnostics(commandArgs(trailingOnly=TRUE),helper_dir)
