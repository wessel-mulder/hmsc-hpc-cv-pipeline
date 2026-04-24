find_model_folders <- function(base_dir = "HmscOutputs", pattern = NULL,
                               full.names = FALSE) {
  folders <- list.dirs(base_dir, recursive = FALSE, full.names = full.names)
  if (!is.null(pattern) && nzchar(pattern)) {
    folders <- folders[grepl(pattern, basename(folders))]
  }
  sort(folders)
}

model_subdir <- function(model_name, ..., base_dir = "HmscOutputs") {
  file.path(base_dir, model_name, ...)
}
