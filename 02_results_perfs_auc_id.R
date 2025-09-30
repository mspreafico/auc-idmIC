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
auc_ID <- weib.auc.id(a01 = a01, a02 = a02, a12 = a12, k = k, t = tab_times)

tab.id.t = NULL
for(lab in c('M','N','O','P','Q','R')){
  
  load(paste0("results/",lab,"_coxtd_id_auc.RData"))
  coxtd_id_t = auc.at.t(object.list = coxtd_id_auc, time.var = "utimes",
                    times = tab_times, true.auc = auc_ID) 
  coxtd_id_t = cbind('Scenario' = lab, 'Model' = 'Cox ROC', coxtd_id_t)
  
  load(paste0("results/",lab,"_mstate_id_auc.RData"))
  excpt_inc = find.except.mstate(mstate_id_auc)
  mstate_id_t = auc.at.t(object.list = mstate_id_auc, except = excpt_inc, inc = TRUE,
                              times = tab_times, true.auc = auc_ID)
  mstate_id_t = cbind('Scenario' = lab, 'Model' = 'Cox prob', mstate_id_t)
  
  load(paste0("results/",lab,"_msm_id_auc.RData"))
  msm_id_t = auc.at.t(object.list = msm_id_auc, times = tab_times, true.auc = auc_ID)
  msm_id_t = cbind('Scenario' = lab, 'Model' = 'PW-const', msm_id_t)
  
  load(paste0("results/",lab,"_smh_id_auc.RData"))
  smh_id_t = auc.at.t(object.list = smh_id_auc, times = tab_times, true.auc = auc_ID)  
  smh_id_t = cbind('Scenario' = lab, 'Model' = 'Weibull', smh_id_t)
  
  load(paste0("results/",lab,"_spline_id_auc.RData"))
  spline_id_t = auc.at.t(object.list = spline_id_auc, times = tab_times, true.auc = auc_ID)  
  spline_id_t = cbind('Scenario' = lab, 'Model' = 'M-spline', spline_id_t)
  
  tab.id.t = rbind.data.frame(tab.id.t, 
                              coxtd_id_t, mstate_id_t, msm_id_t, smh_id_t, spline_id_t)
}

tab.id.t
auc_ID

library(xtable)
tab.id.perfs = cbind.data.frame(tab.id.t[tab.id.t$time==12,c(1:2,6:8)],
                              tab.id.t[tab.id.t$time==36,c(6:8)],
                              tab.id.t[tab.id.t$time==60,c(6:8)])
print(xtable(tab.id.perfs), include.rownames=FALSE)


# make table for incident dynamic validity
tab.id.valid = tab.id.t[tab.id.t$time==60,c(1,2,5)]
tab.id.valid = tab.id.valid[!duplicated(tab.id.valid),]
tab.id.valid
tab.id.valid[tab.id.valid$Model=='M-spline',]

