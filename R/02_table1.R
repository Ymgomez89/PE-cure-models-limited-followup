suppressMessages({library(readxl)})
d <- read_xlsx("data/infarto.xlsx")
egfr <- function(cr, age, sex){k<-ifelse(sex==2,.7,.9);a<-ifelse(sex==2,-.329,-.411);m<-ifelse(sex==2,1.018,1);round(141*pmin((cr/88.4)/k,1)^a*pmax((cr/88.4)/k,1)^(-1.209)*.993^age*m,2)}
d$egfr <- egfr(d$creatinine,d$age,d$sex); g <- d$MACE
sec <- function(s) data.frame(var=s,total="",mace="",non="",p="")
cont <- function(lab,col,dg){ v<-d[[col]]; f<-function(x){q<-quantile(x,c(.5,.25,.75),na.rm=T); fmt<-paste0("%.",dg,"f"); sprintf(paste0(fmt," (",fmt,"–",fmt,")"),q[1],q[2],q[3])}
  p<-wilcox.test(v[g==1],v[g==0],exact=FALSE)$p.value; data.frame(var=lab,total=f(v),mace=f(v[g==1]),non=f(v[g==0]),p=ifelse(p<0.001,"<0.001",sprintf("%.3f",p)))}
catv <- function(lab,col,val){ v<-as.numeric(d[[col]]==val); f<-function(x) sprintf("%d (%.1f)",sum(x),100*mean(x))
  p<-chisq.test(table(v,g))$p.value; data.frame(var=lab,total=f(v),mace=f(v[g==1]),non=f(v[g==0]),p=ifelse(p<0.001,"<0.001",sprintf("%.3f",p)))}
o <- rbind(sec("Clinical characteristics"), cont("Age, years","age",0), catv("Female sex, n (%)","sex",2), cont("Heart rate, beats/min","heart rate",0),
 cont("SBP, mm Hg","systolic blood pressure",0), catv("Anterior wall MI, n (%)","anterior wall myocardial infarction",1), catv("Killip class > I, n (%)","killip grade",2),
 catv("Diabetes mellitus, n (%)","DM",1), catv("Hypertension, n (%)","HTN",1), catv("Previous MI, n (%)","myocardiol infarction history",1),
 sec("Laboratory examinations"), cont("Apelin-12, ng/mL","apelin",2), cont("ΔApelin-12, %","apelin change rate",1), cont("WBC, ×10⁹/L","White blood cells",1),
 cont("Neutrophils, %","neutrophils",1), cont("Haemoglobin, g/L","hemoglobin",0), cont("Albumin, g/L","albumin",1), cont("TC, mmol/L","total cholesterol",2),
 cont("TG, mmol/L","triglycerides",2), cont("HDL-C, mmol/L","High density lipoprotein",2), cont("LDL-C, mmol/L","low density lipoprotein",2),
 cont("FBG, mmol/L","glucose",2), cont("BUN, mmol/L","Urea nitrogen",2), cont("Creatinine, µmol/L","creatinine",1), cont("Uric acid, µmol/L","uric acid",1),
 cont("eGFR, mL/min/1.73 m²","egfr",1), cont("D-dimer, mg/L","D-Dimer",2), cont("Peak CK-MB, U/L","CK-MB",1), cont("Peak cTnI, ng/mL","cTnI",1),
 sec("Cardiac findings"), cont("LAD, mm","left atrial diameter",0), cont("LVEDD, mm","Left ventricular end-diastolic diameter",0),
 catv("Pathological Q wave, n (%)","pathological Q wave",1), cont("GENSINI score","GENSINI score",1))
write.csv(o,"output/Table1.csv",row.names=FALSE); print(o[,c(1,2,5)])
