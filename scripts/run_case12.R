#!/usr/bin/env Rscript
# Entrypoint for Case 1 (banded) + Case 2 (dense)

options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "logs"), showWarnings = FALSE, recursive = TRUE)

set.seed(250)

sink(file.path("results", "logs", "sessionInfo_case12.txt"))
print(sessionInfo())
sink()

source(file.path("scripts", "impl_case12.R"))

cat("Done.\n")
