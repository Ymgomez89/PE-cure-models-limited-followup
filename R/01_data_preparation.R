## 01_data_preparation.R — reads the Dryad dataset and builds the analysis variables
## Input: data/infarto.xlsx (Yang et al. 2017, doi:10.5061/dryad.pf56m)
suppressMessages({library(readxl)})
d0 <- read_xlsx("data/infarto.xlsx")
egfr <- function(cr, age, sex) {           # CKD-EPI 2009; creatinine in µmol/L; sex: 1 = male, 2 = female
  k <- ifelse(sex == 2, 0.7, 0.9); a <- ifelse(sex == 2, -0.329, -0.411); m <- ifelse(sex == 2, 1.018, 1)
  round(141 * pmin((cr / 88.4) / k, 1)^a * pmax((cr / 88.4) / k, 1)^(-1.209) * 0.993^age * m, 2) }
d0$egfr <- egfr(d0$creatinine, d0$age, d0$sex)
## Binary covariates: 1 = present / threshold reached, 0 = otherwise
d <- with(d0, data.frame(time, MACE,
  apelin076 = as.numeric(apelin >= 0.76), killip2 = as.numeric(`killip grade` == 2),
  anterior  = as.numeric(`anterior wall myocardial infarction` == 1), egfr90 = as.numeric(egfr >= 90),
  qwave     = as.numeric(`pathological Q wave` == 1), age63 = as.numeric(age >= 63),
  lad40     = as.numeric(`left atrial diameter` >= 40), lvedd54 = as.numeric(`Left ventricular end-diastolic diameter` >= 54),
  female    = as.numeric(sex == 2), dm = DM, htn = HTN, prevmi = as.numeric(`myocardiol infarction history` == 1)))
f12 <- ~ apelin076 + killip2 + anterior + egfr90 + qwave + age63 + lad40 + lvedd54 + female + dm + htn + prevmi
data12 <- list(y = d$time, delta = d$MACE, z = t(model.matrix(f12, data = d)))
dir.create("output", showWarnings = FALSE)
saveRDS(d0, "output/raw.rds"); saveRDS(d, "output/d12.rds"); saveRDS(data12, "output/data12.rds")
