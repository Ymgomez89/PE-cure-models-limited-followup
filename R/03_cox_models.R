## 03_cox_models.R — Table 2, proportional hazards test, sex-stratified model,
## univariable analyses (Supplementary Table S3) and backward elimination (Supplementary Tables S2 and S4)
suppressMessages(library(survival))
d <- readRDS("output/d12.rds")
vars <- c("apelin076","killip2","anterior","egfr90","qwave","age63","lad40","lvedd54","female","dm","htn","prevmi")
cox <- coxph(as.formula(paste("Surv(time, MACE) ~", paste(vars, collapse = "+"))), data = d)
print(summary(cox)); cat("AIC (partial likelihood):", AIC(cox), "\n")
print(cox.zph(cox))                                                    # scaled Schoenfeld residuals
print(summary(coxph(as.formula(paste("Surv(time, MACE) ~", paste(setdiff(vars,"female"), collapse = "+"), "+ strata(female)")), data = d)))
uni <- t(sapply(vars, function(v) { s <- summary(coxph(as.formula(paste("Surv(time, MACE) ~", v)), data = d)); c(s$conf.int[1, c(1,3,4)], p = s$coefficients[1,5]) }))
write.csv(uni, "output/TableS3_univariable.csv")
cur <- vars                                                            # backward elimination (p > 0.05 removed)
repeat { f <- coxph(as.formula(paste("Surv(time, MACE) ~", paste(cur, collapse = "+"))), data = d)
  p <- summary(f)$coefficients[,5]; if (max(p) <= 0.05) break
  cat("remove", names(which.max(p)), round(max(p),3), "\n"); cur <- setdiff(cur, names(which.max(p))) }
print(summary(f)); cat("Backward-elimination Cox AIC:", AIC(f), "\n")
