#!/usr/bin/env Rscript
# Entrypoint for Case 4 (hub network)

options(stringsAsFactors = FALSE)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("results", "logs"), showWarnings = FALSE, recursive = TRUE)

set.seed(250)

sink(file.path("results", "logs", "sessionInfo_case4_hub.txt"))
print(sessionInfo())
sink()

source(file.path("scripts", "impl_case4_hub.R"))

cat("Done.\n")
