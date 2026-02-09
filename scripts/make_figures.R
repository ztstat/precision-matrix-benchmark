#!/usr/bin/env Rscript
# Figure generation for available results

options(stringsAsFactors = FALSE)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("figures", "plots"), showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
})

save_path <- file.path("figures", "plots")

load_or_stop <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  load(path, envir = .GlobalEnv)
}

get_obj <- function(name) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) get(name, envir = .GlobalEnv) else NULL
}

plot_metric <- function(obj, metric, label_suffix) {
  if (is.null(obj) || is.null(obj[[metric]])) return(invisible(NULL))
  if (is.null(obj$lambda)) return(invisible(NULL))
  
  x <- obj$lambda
  methods <- names(obj[[metric]])
  if (is.null(methods) || length(methods) == 0) return(invisible(NULL))
  
  df <- data.frame(
    x = rep(x, length(methods)),
    y = unlist(obj[[metric]], use.names = FALSE),
    group = rep(methods, each = length(x))
  )
  
  p <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(title = metric, x = expression(lambda), y = metric, color = "Method") +
    theme_minimal()
  
  ggsave(file.path(save_path, paste0(metric, label_suffix, ".png")),
         p, width = 8, height = 6, dpi = 300)
  
  invisible(NULL)
}

plot_bundle <- function(obj, label) {
  for (m in c("BIC", "TPR", "FPR", "Sparsity_ratio", "Fnorm", "Time")) {
    plot_metric(obj, m, label)
  }
}

# Case 1 & 2
if (file.exists(file.path("results", "Results1and2.RData"))) {
  load_or_stop(file.path("results", "Results1and2.RData"))
  for (nm in c("Case1_p100","Case1_p200","Case1_p500","Case2_p100","Case2_p200","Case2_p500")) {
    obj <- get_obj(nm)
    if (!is.null(obj)) plot_bundle(obj, paste0("_", nm))
  }
}

# Case 3
if (file.exists(file.path("results", "Results3_sp0.5.RData"))) {
  load_or_stop(file.path("results", "Results3_sp0.5.RData"))
  for (nm in c("Case3_p100_sp0_5","Case3_p200_sp0_5","Case3_p500_sp0_5")) {
    obj <- get_obj(nm)
    if (!is.null(obj)) plot_bundle(obj, paste0("_", nm))
  }
}
if (file.exists(file.path("results", "Results3_sp0.99.RData"))) {
  load_or_stop(file.path("results", "Results3_sp0.99.RData"))
  for (nm in c("Case3_p100_sp0_99","Case3_p200_sp0_99","Case3_p500_sp0_99")) {
    obj <- get_obj(nm)
    if (!is.null(obj)) plot_bundle(obj, paste0("_", nm))
  }
}

# Case 4
if (file.exists(file.path("results", "Results4_hub0.2.RData"))) {
  load_or_stop(file.path("results", "Results4_hub0.2.RData"))
  for (nm in c("Case4_p100_hub0_2","Case4_p200_hub0_2","Case4_p500_hub0_2")) {
    obj <- get_obj(nm)
    if (!is.null(obj)) plot_bundle(obj, paste0("_", nm))
  }
}

cat("Done.\n")
