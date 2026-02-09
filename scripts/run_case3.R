#!/usr/bin/env Rscript
# Entrypoint for Case 3

options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "logs"), showWarnings = FALSE, recursive = TRUE)

set.seed(250)

sink(file.path("results", "logs", "sessionInfo_case3.txt"))
print(sessionInfo())
sink()

source(file.path("scripts", "impl_case3.R"))

cat("Done.\n")
