#make plots for every 3 scenarios
setwd("~/Library/CloudStorage/OneDrive-UniversiteitLeiden/Leiden/ALICE/home-spreaficom/aucidm")
source('functions/summarize_results.R')
source('functions/weibull_auc.R')

labs = c('M','N','O')
#labs = c('P','Q','R')

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
auc_TID = weib.auc.id(a01 = a01, a02 = a02, a12 = a12, k = k, t = tm)
id_min = 0.45
id_max = round(max(auc_TID, na.rm = TRUE) + 0.1, 1)
true_id = function() lines(tm/12, auc_TID, type = "l", col = "darkblue", lwd = 2)

#######################
# incident dynamic plot
#######################
pdf(paste0("figures/",labs[1], labs[2], labs[3], "_plot_auc_id.pdf"), width = 12, height = 20)
par(mar = c(1, 4.5, 3, 2))
m = matrix(c(1:13, 26, 14:19, 26, 20:26),
            nrow = 7, ncol = 4, byrow = FALSE)
layout(mat = m, 
       widths = c(lcm(5.5), lcm(7.5), lcm(7.5), lcm(7.5)),
       heights = c(lcm(2), lcm(7), lcm(7), lcm(7), 
                   lcm(7), lcm(8.5), lcm(3)))
plot.new()
plot.new()
text(0, 0.5, "Cox ROC", cex = cex.main2, pos = 4)
plot.new()
text(0, 0.5, "Cox prob", cex = cex.main2, pos = 4)
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
load(paste0("results/",lab,"_coxtd_id_auc.RData"))
load(paste0("results/",lab,"_mstate_id_auc.RData"))
excpt_inc = find.except.mstate(mstate_id_auc)
load(paste0("results/",lab,"_msm_id_auc.RData"))
load(paste0("results/",lab,"_smh_id_auc.RData"))
load(paste0("results/",lab,"_spline_id_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "3 months", cex = cex.main2)
plot.all(object.list = coxtd_id_auc, time.var = "utimes",
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = mstate_id_auc,
         title = title, ylim = c(id_min, id_max),
         except = excpt_inc, inc = TRUE)
true_id()
plot.all(object.list = msm_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = smh_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()



# 2. make plots for second scenario
lab = labs[2]
load(paste0("results/",lab,"_coxtd_id_auc.RData"))
load(paste0("results/",lab,"_mstate_id_auc.RData"))
excpt_inc = find.except.mstate(mstate_id_auc)
load(paste0("results/",lab,"_msm_id_auc.RData"))
load(paste0("results/",lab,"_smh_id_auc.RData"))
load(paste0("results/",lab,"_spline_id_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "6 months", cex = cex.main2)
plot.all(object.list = coxtd_id_auc, time.var = "utimes", 
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = mstate_id_auc, 
         title = title, ylim = c(id_min, id_max),
         except = excpt_inc, inc = TRUE)
true_id()
plot.all(object.list = msm_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = smh_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()

# 3. make plots for third scenario
lab = labs[3]
load(paste0("results/",lab,"_coxtd_id_auc.RData"))
load(paste0("results/",lab,"_mstate_id_auc.RData"))
excpt_inc = find.except.mstate(mstate_id_auc)
load(paste0("results/",lab,"_msm_id_auc.RData"))
load(paste0("results/",lab,"_smh_id_auc.RData"))
load(paste0("results/",lab,"_spline_id_auc.RData"))

par(mar = c(1, 4.5, 3, 2))
plot.new()
text(0.5, 0.5, "12 months", cex = cex.main2)
plot.all(object.list = coxtd_id_auc, time.var = "utimes", 
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = mstate_id_auc, 
         title = title, ylim = c(id_min, id_max),
         except = excpt_inc, inc = TRUE)
true_id()
plot.all(object.list = msm_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()
plot.all(object.list = smh_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()

par(mar = c(5, 4.5, 3, 2))
plot.all(object.list = spline_id_auc, 
         title = title, ylim = c(id_min, id_max))
true_id()

# 4. make Legend
par(mar = c(1, 4.5, 1, 2))
plot.new()
legend(x = "top", inset = 0.1,
       legend = c("Truth", "Estimate"), 
       cex = cex.main2, lwd = 6 , col = c("darkblue", rgb(0,1,0,0.3)), 
       horiz = TRUE, bty = "n")

dev.off()
