res <- data.frame(L = 1:20, BerPE = NA, PoPE = NA, NBPE = NA)
for (L in 1:20) for (m in c("BerPE","PoPE","NBPE")) { f <- sprintf("output/scan/%s_%02d.rds", m, L); if (file.exists(f)) res[L, m] <- readRDS(f) }
write.csv(res, "output/AIC_by_partition.csv", row.names = FALSE); print(res)
