# Simulated data for different data scenarios and 
# prepares data for analysis.

source("functions/data_sim_format.R")

n.sim = 1000
N = 400

# scenario M
fup.months = 3
set.seed(20070430)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = FALSE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/M_datasets.RData")


# scenario N
fup.months = 6
set.seed(20070431)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = FALSE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/N_datasets.RData")


# scenario O
fup.months = 12
set.seed(20070432)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = FALSE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/O_datasets.RData")


# scenario P
fup.months = 3
set.seed(20070433)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = TRUE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/P_datasets.RData")


# scenario Q
fup.months = 6
set.seed(20070434)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = TRUE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/Q_datasets.RData")


# scenario R
fup.months = 12
set.seed(20070435)
sim_data <- weib.simu(n_sets = n.sim, n_pt = N, adm.cens = TRUE, ft = fup.months)
coxtd_data <- coxtd.format(sim_data, n_sets = n.sim, n_pt = N)
mstate_data <- mstate.format(sim_data, n_sets = n.sim, n_pt = N)
msm_data <- msm.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)
smh_data <- smh.format(sim_data, n_sets = n.sim, n_pt = N, ft = fup.months)

save(sim_data, coxtd_data, mstate_data, msm_data, smh_data, file = "data/R_datasets.RData")



