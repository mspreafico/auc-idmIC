#make plots for every 3 scenarios
setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/Leiden/ALICE/home-spreaficom/aucidm")
source('functions/summarize_results.R')
source('functions/weibull_auc.R')

#labs = c('M','N','O')
labs = c('P','Q','R')

#######################
# general settings
#######################
# fixed parameters
a01 = 0.05
a02 = 0.05
a12 = 0.56
k = 0.5
w = 60
s = 0:120
t = s + w

# parameters for plotting
lwd1 = 1.5
cex.main1 = 1.5
cex.main2 = 2
title = ""
# compute true AUC lines
tm = seq(0, 120, by = 0.1)
auc_TCD = weib.auc.cd(a01 = a01, a02 = a02, a12 = a12, k = k, s = tm, t = tm + w)
cd_min = round(min(auc_TCD, na.rm = TRUE) - 0.1, 1)
cd_max = round(max(auc_TCD, na.rm = TRUE) + 0.1, 1)
true_cd = function() lines(tm/12, auc_TCD, type = "l", col = "darkblue", lwd = 2)

#######################
# cumulative dynamic plot
#######################
pdf(paste0("figures/",labs[1], labs[2], labs[3], "_plot_auc_cd.pdf"), width = 12, height = 20)
par(mar = c(1, 4.5, 3, 2))
m = matrix(c(1:11, 22, 12:16, 22, 17:22),
           nrow = 6, ncol = 4, byrow = FALSE)
layout(mat = m, 
       widths = c(lcm(5.5), lcm(7.5), lcm(7.5), lcm(7.5)),
       heights = c(lcm(2), lcm(7), lcm(7), lcm(7), lcm(8.5), lcm(3)))
plot.new()
plot.new()
text(0, 0.5, "Cox", cex = cex.main2, pos = 4)
plot.new()
text(0, 0.5, "Piecewise-\nconstant", cex = cex.main2, pos = 4)
plot.new()
text(0, 0.5, "Weibull", cex = cex.main2, pos = 4)

par(mar = c(5, 4.5, 3, 2))
plot.new()
text(0, 0.5, "M-spline", cex = cex.main2, pos = 4)
plot.new()

# 1. make plots for first scenario
lab = labs[1]
load(paste0("results/",lab,"_mstate_cd_auc.RData"))
excpt_cum = find.except.mstate(mstate_cd_auc)
load(paste0("results/",lab,"_msm_cd_auc.RData"))
load(paste0("results/",lab,"_smh_cd_auc.RData"))
load(paste0("results/",lab,"_spline_cd_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "3 months", cex = cex.main2)
plot.all(object.list = mstate_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5),
         except = excpt_cum, inc = FALSE)
true_cd()
plot.all(object.list = msm_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()
plot.all(object.list = smh_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()


# 2. make plots for second scenario
lab = labs[2]
load(paste0("results/",lab,"_mstate_cd_auc.RData"))
excpt_cum = find.except.mstate(mstate_cd_auc)
load(paste0("results/",lab,"_msm_cd_auc.RData"))
load(paste0("results/",lab,"_smh_cd_auc.RData"))
load(paste0("results/",lab,"_spline_cd_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "6 months", cex = cex.main2)
plot.all(object.list = mstate_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5),
         except = excpt_cum, inc = FALSE)
true_cd()
plot.all(object.list = msm_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()
plot.all(object.list = smh_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()

# 3. make plots for third scenario
lab = labs[3]
load(paste0("results/",lab,"_mstate_cd_auc.RData"))
excpt_cum = find.except.mstate(mstate_cd_auc)
load(paste0("results/",lab,"_msm_cd_auc.RData"))
load(paste0("results/",lab,"_smh_cd_auc.RData"))
load(paste0("results/",lab,"_spline_cd_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "12 months", cex = cex.main2)
plot.all(object.list = mstate_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5),
         except = excpt_cum, inc = FALSE)
true_cd()
plot.all(object.list = msm_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()
plot.all(object.list = smh_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_cd_auc, 
         title = title, ylim = c(cd_min, cd_max),
         xlim = c(0, 5))
true_cd()

# 4. make legend
par(mar = c(1, 4.5, 1, 2))
plot.new()
legend(x = "top", inset = 0.1,
       legend = c("Truth", "Estimate"), 
       cex = cex.main2, lwd = 6 , col = c("darkblue", rgb(0,1,0,0.3)), 
       horiz = TRUE, bty = "n")

dev.off()

rm(list=ls())
gc()
