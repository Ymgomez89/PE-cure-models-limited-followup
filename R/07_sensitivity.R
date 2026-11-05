## 07_sensitivity.R — partition sensitivity (Table S5), backward elimination of the PoPE(8) model (Tables S1, S4)
## and nonparametric bootstrap (Table S6; long run time)
suppressMessages(library(survival)); source("R/00_functions.R")
data <- readRDS("output/data12.rds"); d <- readRDS("output/d12.rds")
parts <- list(L4 = part.quantiles(data, 4), L6 = part.quantiles(data, 6), L8 = part.quantiles(data, 8), Clin = c(0, 1, 12, 24))
for (nm in names(parts)) { set.seed(2021); f <- PoPE.est2(data, parts[[nm]])
  cat("\n==", nm, "AIC", f$AIC, "\n"); print(round(f$results, 4)); print(cure.fraction(f, data, nsim = 5000, seed = 2021)) }
vars <- c("apelin076","killip2","anterior","egfr90","qwave","age63","lad40","lvedd54","female","dm","htn","prevmi"); cur <- vars; t8 <- parts$L8
repeat { Z <- t(model.matrix(as.formula(paste("~", paste(cur, collapse = "+"))), data = d))
  set.seed(2021); f <- PoPE.est2(list(y = d$time, delta = d$MACE, z = Z), t8); p <- f$results[2:(length(cur)+1), 5]; names(p) <- cur
  if (max(p) <= 0.05) break; cat("remove", names(which.max(p)), round(max(p), 3), "\n"); cur <- setdiff(cur, names(which.max(p))) }
print(round(f$results, 4)); cat("AIC", f$AIC, "\n")
## Bootstrap: each resample is refitted with exactly the same optimization workflow as the main analysis
## (.fit.optim: Nelder-Mead -> BFGS, convergence code checked at each stage, positive-definite observed
## information), except that the number of restarts is limited to 30 per resample (max.tries = 30).
## Resamples without a converged solution are recorded as failures (NA).
formals(.fit.optim)$max.tries <- 30; n <- length(data$y); boot <- NULL
for (s in 1:4) for (b in 1:50) { set.seed(s * 1000 + b); idx <- sample.int(n, n, replace = TRUE)
  db <- list(y = data$y[idx], delta = data$delta[idx], z = data$z[, idx])
  fb <- tryCatch(PoPE.est2(db, t8), error = function(e) NULL)
  boot <- rbind(boot, if (is.null(fb)) c(rep(NA, 22), NA, NA, NA, NA) else
    c(fb$results[,1], mean(cure.prob(fb, db)), fb$diagnostics$start,
      as.numeric(fb$diagnostics$final.stage == "BFGS"), fb$diagnostics$max.abs.gradient, max(abs(fb$par[1:13])))) }
colnames(boot) <- c(rownames(f0 <- readRDS("output/pe_fits.rds")$PoPE$results), "avg.cure", "starts.used", "final.BFGS", "max.abs.gradient", "max.abs.beta")
write.csv(boot, "output/bootstrap_PoPE8.csv", row.names = FALSE)
cat("converged resamples:", sum(!is.na(boot[, 1])), "of", nrow(boot), "\n")
est <- f0[, 1]; B <- boot[!is.na(boot[, 1]), 1:22]
print(round(cbind(full = c(est, NA), lo = apply(B, 2, quantile, 0.025), hi = apply(B, 2, quantile, 0.975),
  same.sign = c(colMeans(sign(B[, 1:21]) == matrix(sign(est), nrow(B), 21, byrow = TRUE)), NA)), 3))
