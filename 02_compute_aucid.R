on.alice = T
if (on.alice){
  setwd('/home/spreaficom/aucidm')
  cores = 20
}
if (!on.alice) {
  #setwd("~/github/aucidm")
  setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/Leiden/ALICE/home-spreaficom/aucidm")
  cores = parallel::detectCores() - 1
}

library(mcprogress)
lab = 'M'

########################################
# estimate AUC for time-dependent model
########################################
print(paste0(lab,': coxtd model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","coxtd_data")))
source('functions/IncDynAUC.R')

# incident dynamic AUC for time-dependent model
coxtd_id_auc <- pmclapply(coxtd_data, function(x) 
  with(x, risksetAUC.simu(tstop, entry = tstart, status = death, 
                          marker = LR, tmax = max(tstop), plot = FALSE, timefix = FALSE)),
  mc.cores = cores, mc.preschedule = FALSE)
save(coxtd_id_auc, file = paste0("results/",lab,"_coxtd_id_auc.RData"))


########################################
# estimate AUC for mstate model
########################################
print(paste0(lab,': mstate model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","mstate_data")))
load(paste0('models/',lab,'_mstate_model.RData'))
source('functions/IncDynAUC.R')

# incident dynamic AUC 
n = length(mstate_data)
tmat = mstate::transMat(list(c(2, 3), c(3), c()), names = c("ANED", "LR", "Death"))

mstate_id_auc <- pmclapply(1:n, function(j) 
  aucID.mstate(data = mstate_data[[j]], object = mstate_model[[j]], tmat = tmat),
  mc.cores = cores, mc.preschedule = FALSE)
save(mstate_id_auc, file = paste0("results/",lab,"_mstate_id_auc.RData"))


########################################
# estimating AUC for illness-death model
########################################
print(paste0(lab,': msm model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","msm_data")))
load(paste0('models/',lab,'_msm_model.RData'))
source('functions/IncDynAUC.R')

cpoints <- c(6, 30, 60, 90)
lbls <- c(paste("[-Inf, ", cpoints[1], ")", sep=""),
          paste("[", cpoints[1], ", ", cpoints[2],")", sep=""),
          paste("[", cpoints[2], ", ", cpoints[3],")", sep=""),
          paste("[", cpoints[3], ", ", cpoints[4],")", sep=""),
          paste("[", cpoints[4], ", Inf)", sep=""))

# incident dynamic AUC
msm_id_auc <- pmclapply(1:length(msm_data), 
                       function(i) aucID.msm(data = msm_data[[i]], object = msm_model[[i]],
                                             cuts = cpoints, tmp_levels = lbls, marker = "LR"),
                       mc.cores = cores, mc.preschedule = FALSE)
save(msm_id_auc, file = paste0("results/",lab,"_msm_id_auc.RData"))


########################################
# estimating AUC for smooth hazards model with Weibull hazard
########################################
print(paste0(lab,': smh model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","smh_data")))
load(paste0('models/',lab,'_smh_model.RData'))
source('functions/IncDynAUC.R')

# incident dynamic AUC
smh_id_auc <- pmclapply(1:length(smh_data), 
                       function(i) aucID.smh(data = smh_data[[i]], object = smh_model[[i]]),
                       mc.cores = cores, mc.preschedule = FALSE)
save(smh_id_auc, file = paste0("results/",lab,"_smh_id_auc.RData"))


########################################
# estimating AUC for smooth hazard model with splines hazard 
########################################
print(paste0(lab,': spline model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","smh_data")))
load(paste0('models/',lab,'_spline_model.RData'))
source('functions/IncDynAUC.R')

# incident dynamic AUC
spline_id_auc <- pmclapply(1:length(smh_data), 
                          function(i) aucID.smh(data = smh_data[[i]], object = spline_model[[i]]),
                          mc.cores = cores, mc.preschedule = FALSE)
save(spline_id_auc, file = paste0("results/",lab,"_spline_id_auc.RData"))

