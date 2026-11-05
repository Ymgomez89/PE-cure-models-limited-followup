suppressMessages(library(survival)); source("R/00_functions.R")
d <- readRDS("output/d12.rds"); data <- readRDS("output/data12.rds"); out <- readRDS("output/pe_fits.rds"); t8 <- out$PoPE$a
cols <- c("#1b9e77","#d95f02","#7570b3")
dev2 <- function(name, w, h, f) { tiff(paste0("output/",name,".tiff"), width=w, height=h, units="in", res=600, compression="lzw"); f(); dev.off()
  pdf(paste0("output/",name,".pdf"), width=w, height=h); f(); dev.off() }
## Figure 1
f1 <- function(){ km <- survfit(Surv(time, MACE) ~ 1, data=d); par(mar=c(7.5,4.5,1,1))
  plot(km, conf.int=TRUE, ylim=c(0.6,1), xlim=c(0,30), xaxt="n", las=1, xlab="", ylab="MACE-free survival probability", lwd=1.5)
  axis(1, at=seq(0,30,6)); mtext("Time since primary PCI (months)", side=1, line=2.3)
  tt <- seq(0,30,6); nr <- summary(km, times=tt, extend=TRUE)$n.risk
  mtext("Number at risk", side=1, line=4, adj=0, cex=.8, font=2); mtext(nr, side=1, line=5, at=tt, cex=.8) }
dev2("Figure1", 6.5, 4.5, f1)
## Figure 2
f2 <- function(){ x <- read.csv("output/AIC_by_partition.csv"); par(mar=c(4.5,4.5,1,1))
  matplot(x$L, x[,c("BerPE","PoPE","NBPE")], type="b", pch=c(15,16,17), lty=1, col=cols, xlab="Number of intervals (L)", ylab="AIC", las=1, xlim=c(1,20), xaxt="n")
  axis(1, at=1:20); legend("topright", c("BerPE (Bernoulli)","PoPE (Poisson)","NBPE (Negative Binomial)"), pch=c(15,16,17), col=cols, lty=1, bty="n") }
if (file.exists("output/AIC_by_partition.csv")) dev2("Figure2", 6.5, 4.5, f2)
## Figure 3: residual diagnostics (corrected randomized quantile residuals)
po<-out$PoPE; nb<-out$NBPE; be<-out$BerPE; ph<-out$PEPH
cols4 <- c(cols, "grey40")
set.seed(2021)
R <- list(BerPE = nr.BerPE(be$results[,1], data, t8, nrep = 100),
          PoPE  = nr.PoPE(po$results[,1],  data, t8, nrep = 100),
          NBPE  = nr.NBPE(nb$results[,1],  data, t8, nrep = 100),
          "PE-PH (non-cure)" = rqr(Spop.PEPH(ph, data, t8), data$delta, nrep = 100))
rq <- lapply(R, rqr.qq)                                   # median order statistics over 100 randomizations
rs <- t(sapply(R, rqr.summary)); rs <- cbind(rs, max.abs.dev.QQ = sapply(rq, function(r) max(abs(r - qnorm(ppoints(length(r)))))))
print(round(rs, 4)); write.csv(round(rs, 4), "output/TableS10_rqr_summary.csv")
cs <- list(BerPE = -log(Spop.BerPE(be$results[,1], data, t8)), PoPE = -log(Spop.PoPE(po$results[,1], data, t8)),
           NBPE = -log(Spop.NBPE(nb$results[,1], data, t8)), "PE-PH (non-cure)" = -log(Spop.PEPH(ph, data, t8)))
dev.cs <- sapply(cs, function(r) { k <- survfit(Surv(as.numeric(r), data$delta) ~ 1); sel <- k$time <= 0.8
  max(abs(k$cumhaz[sel] - k$time[sel])) }); print(round(dev.cs, 3))
f3 <- function(){ par(mfrow=c(1,2), mar=c(4.5,4.5,2,1))
  qn <- qnorm(ppoints(length(data$y))); rg <- range(c(qn, unlist(rq)))
  plot(NA,xlim=rg,ylim=rg,xlab="N(0,1) quantiles",ylab="Randomized quantile residuals",main="(a) Q-Q plot",las=1)
  for(i in 1:4) points(qn, rq[[i]], col=cols4[i], pch=c(16,16,16,1)[i], cex=.5); abline(0,1,lty=2)
  legend("topleft",names(rq),col=cols4,pch=c(16,16,16,1),bty="n",cex=.8)
  plot(NA,xlim=c(0,1.15),ylim=c(0,1.3),xlab="Cox-Snell residual",ylab="Estimated cumulative hazard",main="(b) Cox-Snell residuals",las=1)
  for(i in 1:4){k<-survfit(Surv(as.numeric(cs[[i]]),data$delta)~1); lines(k$time,k$cumhaz,type="s",col=cols4[i],lwd=1.5,lty=c(1,1,1,2)[i])}
  abline(0,1,lty=3); legend("topleft",names(cs),col=cols4,lwd=1.5,lty=c(1,1,1,2),bty="n",cex=.8) }
dev2("Figure3", 9, 4.5, f3)
## Supplementary Figure S1: profile likelihood
prof <- read.csv("output/profile_cure_fraction.csv")
fs1 <- function(){ par(mar=c(4.5,4.5,1,1)); plot(prof$cure, prof$LRT, type="l", lwd=2, xlab="Average cure fraction", ylab="Likelihood-ratio statistic", las=1, ylim=c(0,16))
  abline(h=qchisq(.95,1), lty=2); abline(v=0.66, lty=3); text(0.12, qchisq(.95,1)+0.6, "95% threshold (3.84)", adj=0, cex=.8) }
dev2("FigureS1", 6, 4.5, fs1)
png("output/FigureS1.png", width=1800, height=1350, res=300); fs1(); dev.off()
