
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
w.hor = 60
s.vec = 0:120

########################################
# estimate AUC for mstate model
########################################
print(paste0(lab,': mstate model'))
load(paste0('data/',lab,'_datasets.RData'))
rm(list=setdiff(ls(), c("lab","cores","w.hor","s.vec","mstate_data")))
load(paste0('models/',lab,'_mstate_model.RData'))
source('functions/CumDynAUC.R')

# cumulative dynamic AUC 
n = length(mstate_data)
tmat = mstate::transMat(list(c(2, 3), c(3), c()), names = c("ANED", "LR", "Death"))
mstate_cd_auc <- pmclapply(1:n, function(j) 
  aucCD.mstate(data = mstate_data[[j]], object = mstate_model[[j]], 
               tmat = tmat, s = s.vec, w = w.hor),
  mc.cores = cores, mc.preschedule = FALSE)
save(mstate_cd_auc, file = paste0("results/",lab,"_mstate_cd_auc.RData"))


########################################
# estimating AUC for illness-death model
########################################
print(paste0(lab,': msm model'))
rm(list=setdiff(ls(), c("lab","cores","w.hor","s.vec")))
load(paste0('models/',lab,'_msm_model.RData'))
source('functions/CumDynAUC.R')

cpoints <- c(6, 30, 60, 90)
lbls <- c(paste("[-Inf, ", cpoints[1], ")", sep=""),
          paste("[", cpoints[1], ", ", cpoints[2],")", sep=""),
          paste("[", cpoints[2], ", ", cpoints[3],")", sep=""),
          paste("[", cpoints[3], ", ", cpoints[4],")", sep=""),
          paste("[", cpoints[4], ", Inf)", sep=""))

# cumulative dynamic AUC
# data can be in wide or long format
msm_cd_auc <- pmclapply(msm_model, 
                        function(x) aucCD.msm(object = x, s=s.vec, w=w.hor,
                                              tmp_levels = lbls, marker = "LR"),
                        mc.cores = cores, mc.preschedule = FALSE)
save(msm_cd_auc, file = paste0("results/",lab,"_msm_cd_auc.RData"))


########################################
# estimating AUC for smooth hazards model with weibul hazard
########################################
print(paste0(lab,': smh model'))
rm(list=setdiff(ls(), c("lab","cores","w.hor","s.vec")))
load(paste0('models/',lab,'_smh_model.RData'))
source('functions/CumDynAUC.R')

# cumulative dynamic AUC
smh_cd_auc <- pmclapply(smh_model, 
                       function(x) aucCD.smh(object = x, s=s.vec, w=w.hor),
                       mc.cores = cores, mc.preschedule = FALSE)
save(smh_cd_auc, file = paste0("results/",lab,"_smh_cd_auc.RData"))


########################################
# estimating AUC for smooth hazard model with spline hazard 
########################################
print(paste0(lab,': spline model'))
rm(list=setdiff(ls(), c("lab","cores","w.hor","s.vec")))
load(paste0('models/',lab,'_spline_model.RData'))
source('functions/CumDynAUC.R')

# cumulative dynamic AUC
spline_cd_auc <- pmclapply(spline_model, 
                          function(x) aucCD.smh(object = x, s = s.vec, w = w.hor),
                          mc.cores = cores, mc.preschedule = FALSE)
save(spline_cd_auc, file = paste0("results/",lab,"_spline_cd_auc.RData"))
