# Utilities: seeding, timing, saving outputs

with_seed <- function(seed, expr) {
  if (!is.null(seed)) set.seed(seed)
  force(expr)
}

save_results_rds <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
}
