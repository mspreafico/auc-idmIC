setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/Leiden/ALICE/home-spreaficom/aucidm")

source('functions/summarize_results.R')

# Weibull parameters
a01 = 0.05
a02 = 0.05
a12 = 0.56
k = 0.5
true_coef = log(a12/a02)

tab = NULL
for(lab in c('M','N','O','P','Q','R')){
  
  load(paste0("models/",lab,"_coxtd_model.RData"))
  cox_perfs = model.perf(coxtd_model, true.coef = true_coef, msm.type=F)
  load(paste0("models/",lab,"_mstate_model.RData"))
  mstate_perfs = model.perf(mstate_model, true.coef = true_coef, msm.type=F)
  load(paste0("models/",lab,"_msm_model.RData"))
  msm_perfs = model.perf(msm_model, true.coef = true_coef, msm.type=T)
  
  tab = rbind.data.frame(tab,
                         cbind.data.frame(
                           'Scenario' = rep(lab,3),
                           'Model' = c('Cox ROC','Cox prob','PW-const'),
                           'coef' = c(cox_perfs$LR,mstate_perfs$LR,msm_perfs$LR),
                           'exp(coef)' = c(cox_perfs$LR_HR,mstate_perfs$LR_HR,msm_perfs$LR_HR),
                           'empSE' = c(cox_perfs$seLR,mstate_perfs$seLR,msm_perfs$seLR),
                           'Bias' = c(cox_perfs$bias,mstate_perfs$bias,msm_perfs$bias),
                           'RMSE' = c(cox_perfs$rmse,mstate_perfs$rmse,msm_perfs$rmse)
                         ))
  
}
tab

save(tab, file = "results/results_MR_coefs.RData")

library(xtable)
print(xtable(tab), include.rownames=FALSE)

