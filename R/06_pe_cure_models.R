## 06_pe_cure_models.R — PE cure models with L = 8 (Table 4), cure fraction, profile likelihood (Figure S1),
## events / patients at risk / person-time per interval (Supplementary Table S7), non-cure PE proportional
## hazards comparator with convergence diagnostics and Poisson-GLM verification (Table 5; Supplementary Table S8),
## multiple-start verification of all PE models (Supplementary Table S9)
suppressMessages(library(survival)); source("R/00_functions.R")
data <- readRDS("output/data12.rds"); d <- readRDS("output/d12.rds"); t8 <- part.quantiles(data, 8); fits <- list()
s <- nrow(data$z); vars <- rownames(data$z)[-1]
for (m in c("PoPE","NBPE","BerPE")) { set.seed(2021); fits[[m]] <- get(paste0(m,".est2"))(data, t8) }
## non-cure PE proportional hazards comparator (same covariates and intervals)
set.seed(2021); fits$PEPH <- PEPH.est2(data, t8)
po <- fits$PoPE; saveRDS(fits, "output/pe_fits.rds")
print(round(po$results, 4)); cat("logLik", po$log.likelihood, "k", length(po$par), "AIC", po$AIC, "\n")
print(cure.fraction(po, data, nsim = 5000, seed = 2021)); print(range(cure.prob(po, data)))
## Supplementary Table S7: events, patients at risk at the start of each interval and person-time
ev <- events.by.interval(data, t8); brk <- c(t8, Inf)
ev$person.time <- sapply(1:8, function(l) sum(pmax(0, pmin(data$y, brk[l + 1]) - brk[l])))
ev$crude.rate <- ev$events / ev$person.time
write.csv(ev, "output/TableS7_events_by_interval.csv", row.names = FALSE); print(ev)
cat("events per estimated parameter, PoPE(8):", sum(data$delta) / length(po$par), "\n")
## convergence diagnostics of the final fits
diagtab <- do.call(rbind, lapply(names(fits), function(m) data.frame(model = m, k = length(fits[[m]]$par),
  logLik = fits[[m]]$log.likelihood, AIC = fits[[m]]$AIC, as.data.frame(fits[[m]]$diagnostics))))
print(diagtab); write.csv(diagtab, "output/convergence_diagnostics.csv", row.names = FALSE)
## comparator: estimates, hazard ratios and Poisson-GLM verification
pe <- fits$PEPH; print(round(pe$results, 5)); print(round(pe$HR, 4))
gl <- PEPH.glm(d, vars, t8)
chk <- cbind(gamma.ML = pe$par[1:(s - 1)], gamma.GLM = coef(gl$fit)[vars],
             se.ML = sqrt(diag(pe$vcov))[1:(s - 1)], se.GLM = sqrt(diag(vcov(gl$fit)))[vars])
print(round(chk, 6)); cat("max |diff| coef:", max(abs(chk[,1] - chk[,2])), " max |diff| SE:", max(abs(chk[,3] - chk[,4])), "\n")
cat("PE-PH logLik (optim):", pe$log.likelihood, "  (Poisson GLM):", gl$logLik.PEPH, "\n")
cat("max |diff| log lambda:", max(abs(pe$par[s:(s + 7)] - coef(gl$fit)[1:8])), "\n")
write.csv(data.frame(round(pe$HR, 4), round(chk, 6)), "output/TableS8_PEPH_estimates.csv")
write.csv(round(pe$results, 6), "output/TableS8_PEPH_full_results.csv")
## multiple-start verification (50 random starts per model, perturbation SD 0.5 on the internal scale)
ini0 <- function(extra = 0) function() c(rep(0, s), rep(-log(mean(data$y)), 8), rep(0, extra)) + rnorm(s + 8 + extra, sd = 0.5)
iniPH <- function() c(rep(0, s - 1), rep(log(sum(data$delta) / sum(data$y)), 8)) + rnorm(s - 1 + 8, sd = 0.5)
ms <- list(PoPE = multistart.check(llike.PoPE, ini0(), data, t8),
           NBPE = multistart.check(llike.NBPE, ini0(1), data, t8),
           BerPE = multistart.check(llike.BerPE, ini0(), data, t8),
           PEPH = multistart.check(llike.PEPH, iniPH, data, t8))
mstab <- do.call(rbind, lapply(names(ms), function(m) { x <- ms[[m]]; best <- fits[[m]]$log.likelihood
  data.frame(model = m, starts = nrow(x), converged = sum(!is.na(x$logLik)),
    best.logLik = max(x$logLik, na.rm = TRUE), reported.logLik = best,
    within.0.001 = sum(abs(x$logLik - max(x$logLik, na.rm = TRUE)) < 1e-3, na.rm = TRUE),
    min.logLik = min(x$logLik, na.rm = TRUE), max.grad = max(x$max.grad, na.rm = TRUE)) }))
print(mstab); write.csv(mstab, "output/TableS9_multistart.csv", row.names = FALSE); saveRDS(ms, "output/multistart.rds")
## profile likelihood of the average cure fraction
Z <- data$z
avgcure <- function(b) mean(plogis(drop(t(Z) %*% b)))
b0_for <- function(c, br) uniroot(function(b0) avgcure(c(b0, br)) - c, c(-30, 30))$root
prof <- NULL; set.seed(2021)
for (c in seq(0.10, 0.80, by = 0.025)) {
  nll <- function(p, data, t) { b0 <- tryCatch(b0_for(c, p[1:(s-1)]), error = function(e) NA); if (is.na(b0)) return(1e10); llike.PoPE(c(b0, p), data, t) }
  best <- NULL
  for (k in 1:3) { ini <- po$par[-1] + rnorm(length(po$par) - 1, sd = 0.05)
    r <- try(optim(ini, nll, data = data, t = t8, method = "BFGS", control = list(maxit = 3000)), silent = TRUE)
    if (!inherits(r, "try-error") && r$convergence == 0 && (is.null(best) || r$value < best$value)) best <- r }
  prof <- rbind(prof, data.frame(cure = c, LRT = 2 * (best$value + po$log.likelihood))) }
write.csv(prof, "output/profile_cure_fraction.csv", row.names = FALSE); print(prof)
