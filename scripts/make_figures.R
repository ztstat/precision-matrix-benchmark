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

plot_metric <- function(obj, metric, label_suffix, x) {
  if (is.null(obj) || is.null(obj[[metric]])) return(invisible(NULL))
  
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

plot_bundle <- function(obj, label, x) {
  for (m in c("BIC", "TPR", "FPR", "Sparsity_ratio", "Fnorm", "Time")) {
    plot_metric(obj, m, label, x = x)
  }
}

detect_xgrid <- function(obj) {
  # Default quick grid
  x <- seq(0.2, 1.5, length.out = 10)
  
  # If any metric vector length differs, adapt x to match
  pick <- NULL
  for (m in c("BIC","TPR","FPR","Sparsity_ratio","Fnorm","Time")) {
    if (!is.null(obj[[m]]) && length(obj[[m]]) > 0) {
      v <- obj[[m]][[1]]
      if (!is.null(v) && length(v) > 0) {
        pick <- length(v)
        break
      }
    }
  }
  if (!is.null(pick) && pick != length(x)) {
    x <- seq(0.1, 2, length.out = pick)
  }
  x
}

# -----------------------------
# Case 1 & 2
# -----------------------------
if (file.exists(file.path("results", "Results1and2.RData"))) {
  load_or_stop(file.path("results", "Results1and2.RData"))
  
  candidates <- c(
    "Case1_p100","Case1_p200","Case1_p500",
    "Case2_p100","Case2_p200","Case2_p500"
  )
  
  for (nm in candidates) {
    obj <- get_obj(nm)
    if (is.null(obj)) next
    x <- detect_xgrid(obj)
    plot_bundle(obj, label = paste0("_", nm), x = x)
  }
}

# -----------------------------
# Case 3
# -----------------------------
if (file.exists(file.path("results", "Results3_0.5.RData"))) {
  load_or_stop(file.path("results", "Results3_0.5.RData"))
}
if (file.exists(file.path("results", "Results3_0.99.RData"))) {
  load_or_stop(file.path("results", "Results3_0.99.RData"))
}

case3_candidates <- c(
  "Case3_p=100_Sp=0.5",
  "Case3_p=200_Sp=0.5",
  "Case3_p=500_Sp=0.5",
  "Case3_p=100_Sp=0.99",
  "Case3_p=200_Sp=0.99",
  "Case3_p=500_Sp=0.99"
)

for (nm in case3_candidates) {
  obj <- get_obj(nm)
  if (is.null(obj)) next
  x <- detect_xgrid(obj)
  plot_bundle(obj, label = paste0("_", gsub("[^A-Za-z0-9_=-]", "", nm)), x = x)
}

# -----------------------------
# Case 4
# -----------------------------
if (file.exists(file.path("results", "Results_4_hub=0.2.RData"))) {
  load_or_stop(file.path("results", "Results_4_hub=0.2.RData"))
  
  case4_candidates <- c(
    "Case4_p=100_hub=0.2",
    "Case4_p=200_hub=0.2",
    "Case4_p=500_hub=0.2"
  )
  
  for (nm in case4_candidates) {
    obj <- get_obj(nm)
    if (is.null(obj)) next
    x <- detect_xgrid(obj)
    plot_bundle(obj, label = paste0("_", gsub("[^A-Za-z0-9_=-]", "", nm)), x = x)
  }
}

cat("Done.\n")
