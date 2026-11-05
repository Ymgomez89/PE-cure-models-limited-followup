## Optional helper: runs one cell (L, model) of 05_pe_partition_scan.R with identical settings.
## Usage: Rscript R/05b_partition_scan_chunk.R L model ; results in output/scan/; assembled by 05c.
suppressMessages(library(survival)); source("R/00_functions.R")
a <- commandArgs(TRUE); L <- as.integer(a[1]); m <- a[2]
formals(.fit.optim)$max.tries <- 30   # at most 30 starts per fit in the partition scan
data <- readRDS("output/data12.rds"); t <- part.quantiles(data, L)
set.seed(2021); aic <- tryCatch(get(paste0(m, ".est2"))(data, t)$AIC, error = function(e) NA)
saveRDS(aic, sprintf("output/scan/%s_%02d.rds", m, L))
