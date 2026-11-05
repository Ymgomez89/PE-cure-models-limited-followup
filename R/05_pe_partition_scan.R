## 05_pe_partition_scan.R — AIC of BerPE, PoPE and NBPE for L = 1,...,20 (Figure 2). Long run time.
suppressMessages(library(survival)); source("R/00_functions.R")
formals(.fit.optim)$max.tries <- 30   # at most 30 starts per fit in the partition scan
data <- readRDS("output/data12.rds"); res <- data.frame(L = 1:20, BerPE = NA, PoPE = NA, NBPE = NA)
for (L in 1:20) { t <- part.quantiles(data, L)
  for (m in c("BerPE","PoPE","NBPE")) { set.seed(2021); res[L, m] <- tryCatch(get(paste0(m,".est2"))(data, t)$AIC, error = function(e) NA) } }
write.csv(res, "output/AIC_by_partition.csv", row.names = FALSE)
