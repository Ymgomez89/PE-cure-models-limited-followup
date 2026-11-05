## Helper used to run the bootstrap of 07_sensitivity.R in chunks (identical seeds and settings).
## Usage: Rscript R/07b_bootstrap_chunk.R s b_from b_to
suppressMessages(library(survival)); source("R/00_functions.R")
a <- as.integer(commandArgs(TRUE)); s <- a[1]
data <- readRDS("output/data12.rds"); t8 <- part.quantiles(data, 8)
formals(.fit.optim)$max.tries <- 30; n <- length(data$y)
for (b in a[2]:a[3]) { f <- sprintf("output/boot/b_%d_%02d.rds", s, b); if (file.exists(f)) next
  set.seed(s * 1000 + b); idx <- sample.int(n, n, replace = TRUE)
  db <- list(y = data$y[idx], delta = data$delta[idx], z = data$z[, idx])
  fb <- tryCatch(PoPE.est2(db, t8), error = function(e) NULL)
  saveRDS(if (is.null(fb)) c(rep(NA, 22), NA, NA, NA, NA) else
    c(fb$results[,1], mean(cure.prob(fb, db)), fb$diagnostics$start,
      as.numeric(fb$diagnostics$final.stage == "BFGS"), fb$diagnostics$max.abs.gradient, max(abs(fb$par[1:13]))), f) }
