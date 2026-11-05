## 04_weibull_cure.R — Weibull cure models (Table 3)
## Models whose likelihood is maximized at the boundary of the parameter space (diverging coefficients or
## scale parameter; no finite MLE with a positive-definite observed information matrix) are reported as
## not identifiable (n.i.): their supremum of the log-likelihood is obtained from 20 random starts
## (Nelder-Mead followed by BFGS; the convergence codes of the best solution are reported).
suppressMessages(library(survival)); source("R/00_functions.R")
data <- readRDS("output/data12.rds"); s <- nrow(data$z); out <- list()
sup.fit <- function(nll, ini.fun, nstart = 20) { best <- NULL
  for (k in 1:nstart) { r1 <- try(optim(ini.fun(), nll, data = data, t = 0, control = list(maxit = 10000)), silent = TRUE)
    if (inherits(r1, "try-error") || !is.finite(r1$value)) next
    r1$codes <- c(NM = r1$convergence, BFGS = NA)
    r2 <- try(optim(r1$par, nll, data = data, t = 0, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(r2, "try-error")) { r1$codes["BFGS"] <- r2$convergence
      if (is.finite(r2$value) && r2$value <= r1$value) { cd <- r1$codes; r1 <- r2; r1$codes <- cd } }
    if (is.null(best) || r1$value < best$value) best <- r1 }
  best }
formals(.fit.optim)$max.tries <- 30   # non-identifiable models fail all starts; 30 starts suffice here
for (m in c("BerW","PoW","NBW")) { set.seed(2021)
  f <- tryCatch(get(paste0(m,".est"))(data), error = function(e) NULL)
  if (is.null(f)) { nll <- get(paste0("llike.", m)); k <- s + 2 + (m == "NBW")
    ini <- if (m == "NBW") function() c(.ini.w(data, s), log(1.5)) + c(rep(0, s + 2), rnorm(1, sd = 0.1)) else function() .ini.w(data, s)
    set.seed(2021); b <- sup.fit(nll, ini)
    f <- list(log.likelihood = -b$value, AIC = 2 * b$value + 2 * k, par = b$par, identifiable = FALSE,
              alpha = exp(b$par[s + 1]), lambda = exp(b$par[s + 2]), max.abs.beta = max(abs(b$par[1:s])), codes = b$codes) } else
    f$identifiable <- max(abs(f$par[1:s])) < 10 && f$results["lambda", 1] < 1e4
  out[[m]] <- f
  cat("\n==", m, " logLik", out[[m]]$log.likelihood, " AIC", out[[m]]$AIC, " identifiable:", out[[m]]$identifiable, "\n")
  if (!is.null(f$results)) print(round(f$results, 4)) else print(c(alpha = f$alpha, lambda = f$lambda, max.abs.beta = f$max.abs.beta, f$codes)) }
saveRDS(out, "output/weibull_fits.rds")
