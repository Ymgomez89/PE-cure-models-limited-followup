## 07c_bootstrap_summary.R — Supplementary Table S6. Uses output/bootstrap_PoPE8.csv (from 07_sensitivity.R)
## or, if absent, assembles the chunk files written by 07b_bootstrap_chunk.R.
suppressMessages(library(survival)); source("R/00_functions.R")
f0 <- readRDS("output/pe_fits.rds")$PoPE; data <- readRDS("output/data12.rds")
nm <- c(rownames(f0$results), "avg.cure", "starts.used", "final.BFGS", "max.abs.gradient", "max.abs.beta")
if (file.exists("output/bootstrap_PoPE8.csv")) boot <- as.matrix(read.csv("output/bootstrap_PoPE8.csv")) else {
  boot <- do.call(rbind, lapply(1:4, function(s) do.call(rbind, lapply(1:50, function(b) readRDS(sprintf("output/boot/b_%d_%02d.rds", s, b))))))
  colnames(boot) <- nm; write.csv(boot, "output/bootstrap_PoPE8.csv", row.names = FALSE) }
ok <- !is.na(boot[, 1]); cat("converged resamples:", sum(ok), "of", nrow(boot), "\n")
cat("final solution from BFGS:", sum(boot[ok, "final.BFGS"] == 1), "; starts used: median", median(boot[ok, "starts.used"]),
    "max", max(boot[ok, "starts.used"]), "; max |gradient|:", signif(max(boot[ok, "max.abs.gradient"]), 3), "\n")
B <- boot[ok, 1:22]; est <- c(f0$results[, 1], mean(cure.prob(f0, data)))
tab <- cbind(full = est, median = apply(B, 2, median), lo = apply(B, 2, quantile, 0.025), hi = apply(B, 2, quantile, 0.975),
  same.sign = c(100 * colMeans(sign(B[, 1:13]) == matrix(sign(est[1:13]), nrow(B), 13, byrow = TRUE)), rep(NA, 9)))
print(round(tab, 3)); write.csv(round(tab, 4), "output/TableS6_bootstrap_summary.csv")
## resamples with divergent cure-component coefficients (monotone likelihood; no finite MLE)
dv <- boot[ok, "max.abs.beta"] > 10
cat("divergent resamples (max |beta| > 10):", sum(dv), "\n")
Bf <- B[!dv, ]
tabf <- cbind(median = apply(Bf, 2, median), lo = apply(Bf, 2, quantile, 0.025), hi = apply(Bf, 2, quantile, 0.975),
  same.sign = c(100 * colMeans(sign(Bf[, 1:13]) == matrix(sign(est[1:13]), nrow(Bf), 13, byrow = TRUE)), rep(NA, 9)))
cat("Excluding divergent resamples:\n"); print(round(tabf, 3)); write.csv(round(tabf, 4), "output/TableS6_bootstrap_summary_finite.csv")
