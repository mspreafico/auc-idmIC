setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/Leiden/ALICE/home-spreaficom/aucidm")

source('functions/summarize_results.R')
source('functions/weibull_auc.R')

# Weibull parameters
a01 = 0.05
a02 = 0.05
a12 = 0.56
k = 0.5

#######################
# compute average  auc at specific time points (tab_times)
#######################
w.hor = 60
tab_times = c(1, 3, 5)*12

# Incident Dynamic AUC
auc_CD <- weib.auc.cd(a01 = a01, a02 = a02, a12 = a12, k = k, 
                      s = tab_times, t = tab_times + w.hor)

tab.cd.t = NULL
for(lab in c('M','N','O','P','Q','R')){
  
  load(paste0("results/",lab,"_mstate_cd_auc.RData"))
  excpt_cum = find.except.mstate(mstate_cd_auc)
  mstate_cd_t = auc.at.t(object.list = mstate_cd_auc, except = excpt_cum, inc = FALSE,
                         times = tab_times, true.auc = auc_CD)
  mstate_cd_t = cbind('Scenario' = lab, 'Model' = 'Cox prob', mstate_cd_t)
  
  load(paste0("results/",lab,"_msm_cd_auc.RData"))
  msm_cd_t = auc.at.t(object.list = msm_cd_auc, times = tab_times, true.auc = auc_CD)
  msm_cd_t = cbind('Scenario' = lab, 'Model' = 'PW-const', msm_cd_t)
  
  load(paste0("results/",lab,"_smh_cd_auc.RData"))
  smh_cd_t = auc.at.t(object.list = smh_cd_auc, times = tab_times, true.auc = auc_CD)  
  smh_cd_t = cbind('Scenario' = lab, 'Model' = 'Weibull', smh_cd_t)
  
  load(paste0("results/",lab,"_spline_cd_auc.RData"))
  spline_cd_t = auc.at.t(object.list = spline_cd_auc, times = tab_times, true.auc = auc_CD)  
  spline_cd_t = cbind('Scenario' = lab, 'Model' = 'M-spline', spline_cd_t)
  
  tab.cd.t = rbind.data.frame(tab.cd.t, 
                              mstate_cd_t, msm_cd_t, smh_cd_t, spline_cd_t)
}
tab.cd.t
auc_CD

library(xtable)
tab.cd.perfs = cbind.data.frame(tab.cd.t[tab.cd.t$time==12,c(1:2,6:8)],
                             tab.cd.t[tab.cd.t$time==36,c(6:8)],
                             tab.cd.t[tab.cd.t$time==60,c(6:8)])
print(xtable(tab.cd.perfs), include.rownames=FALSE)

# make table for cumulative dynamic validity
tab.cd.valid = tab.cd.t[tab.cd.t$time==60,c(1,2,5)]
tab.cd.valid = tab.cd.valid[!duplicated(tab.cd.valid),]
tab.cd.valid
tab.cd.valid[tab.cd.valid$Model=='M-spline',]
