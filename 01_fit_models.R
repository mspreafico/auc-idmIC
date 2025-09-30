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

#############################
# fit models 
#############################

# fit time-dependent Cox model
library(survival)
print(paste0(lab,': coxtd model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","coxtd_data")))

frml <- Surv(tstart, tstop, death) ~ LR
coxtd_model <- pmclapply(coxtd_data, function(x) coxph(frml, data = x, timefix = FALSE), 
                        mc.cores = cores, mc.preschedule = FALSE) 
save(coxtd_model, file = paste0("models/",lab,"_coxtd_model.RData"))


# fit mstate model
library(mstate)
print(paste0(lab,': mstate model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","mstate_data")))

frml <- Surv(Tstart, Tstop, status) ~ disease + strata(strata)
mstate_model <- pmclapply(mstate_data, function(x) coxph(frml, data = x, timefix = FALSE), 
                         mc.cores = cores, mc.preschedule = FALSE) 
save(mstate_model, file = paste0("models/",lab,"_mstate_model.RData"))


# fit interval censored illness-death model with msm
library(msm)
print(paste0(lab,': msm model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","msm_data")))
source("functions/fit_msm.R")

qmat <- rbind(c(0, 0.3, 0.3), c(0, 0, 0.3), 
              c(0, 0, 0)) 
rownames(qmat) <- colnames(qmat) <- c("ANED", "LR", "Death")
#msm_model <- mclapply(msm_data, function(x) fit_msm(x), mc.cores = cores, mc.preschedule = FALSE)
msm_model <- pmclapply(msm_data, function(x) fit_msm(x), 
                       mc.cores = cores, mc.preschedule = FALSE)
save(msm_model, file = paste0("models/",lab,"_msm_model.RData"))


#  fit smooth hazards model with weibul hazard
library(SmoothHazard)
print(paste0(lab,': smh model'))
load(paste0("data/",lab,"_datasets.RData"))
rm(list=setdiff(ls(), c("lab","cores","smh_data")))

smh_model <- pmclapply(smh_data, 
                       function(x) idm(formula01 = Hist(time = list(L, R), event = d1) ~ 1, 
                                       formula02 = Hist(time = tt, event = d2) ~ 1, 
                                       formula12 = ~ 1, data = x, method="Weib"), 
                       mc.cores = cores, mc.preschedule = FALSE)
save(smh_model, file = paste0("models/",lab,"_smh_model.RData"))


# fit smooth hazards model with spline hazard
print(paste0(lab,': spline model'))
spline_model <- pmclapply(smh_data, 
                          function(x) idm(formula01 = Hist(time = list(L, R), event = d1) ~ 1, 
                                          formula02 = Hist(time = tt, event = d2) ~ 1, 
                                          formula12 = ~ 1, data = x, method = "Splines"), 
                          mc.cores = cores, mc.preschedule = FALSE)
save(spline_model, file = paste0("models/",lab,"_spline_model.RData"))


