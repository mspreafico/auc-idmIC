library(mstate)
library(msm)
library(SmoothHazard)
library(risksetROC)

##############################
# modification of H&Z function to include timefix option of coxph
# (needed because data is simulated)
##############################

risksetAUC.simu <- function (Stime, entry = NULL, status, marker, method = "Cox", 
                        span = NULL, order = 1, window = "asymmetric", tmax, weight = "rescale", 
                        plot = TRUE, type = "l", xlab = "Time", ylab = "AUC", timefix=FALSE, ...) 
{
  mChoice <- match(method, c("Cox", "LocalCox", "Schoenfeld"))
  if (is.na(mChoice)) {
    cat("error in method choice")
    stop(0)
  }
  if (is.null(span) & ((mChoice == 2) || (mChoice == 3))) {
    cat("Need span for methods = \"LocalCox\" or \"Schoenfeld\" \n")
    stop(0)
  }
  p = 1
  eta = marker
  if (length(entry) == 0) {
    time = rep(0, length(Stime))
  }
  else {
    time = entry
  }
  Stime = Stime
  time2 = Stime
  utimes = unique(Stime[status == 1])
  utimes = utimes[order(utimes)]
  new.eta = NULL
  km.out = weightedKM(Stime = Stime, status = status, entry = time)
  if (mChoice == 1) {
    fit = coxph(Surv(time = time, time2 = time2, event = status) ~ 
                  eta, timefix=timefix)
    gamma.t = rep(fit$coefficients, NROW(utimes))
  }
  if (mChoice == 2) {
    bfnx.ll = llCoxReg(Stime = Stime, status = status, marker = eta, 
                       span = span, p = p, window = window, entry = entry)
    gamma.t = bfnx.ll$beta[, 1]
  }
  if (mChoice == 3) {
    fit = coxph(Surv(time = time, time2 = time2, event = status) ~ 
                  eta, timefix=timefix)
    bfnx.SS = SchoenSmooth(fit = fit, Stime = Stime, status = status, 
                           span = span, order = order)
    gamma.t = bfnx.SS$beta
    gamma.t = gamma.t[!duplicated(gamma.t)]
  }
  AUC = NULL
  for (i in 1:NROW(gamma.t)) {
    new.eta = eta * gamma.t[i]
    out = CoxWeights(marker = new.eta, Stime = Stime, status = status, 
                     predict.time = utimes[i], entry = entry)
    AUC = c(AUC, out$AUC)
  }
  Cindex = IntegrateAUC(AUC = AUC, utimes = utimes, St = km.out$survival, 
                        tmax = tmax, weight = weight)
  if (plot == TRUE) {
    plot(utimes, AUC, type = type, xlim = c(min(utimes), 
                                            tmax + 1), ylim = c(0.4, 1), xlab = xlab, ylab = ylab, 
         ...)
    abline(h = 0.5)
  }
  return(out = list(utimes = utimes, St = km.out$survival, 
                    AUC = AUC, Cindex = Cindex))
}

##################################
# functions to estimate AUC_ID(t)
##################################
# estimates incidence dynamic AUC for mstate model
# x is data, needed for msfit
aucID.mstate <- function(data, object, tmat){
    assign("x", data, envir = .GlobalEnv)
    ms <- msfit(object, newdata = data.frame(disease = c(0, 0, 1),
                                             trans = 1:3,
                                             strata = c(1, 2, 2)),
                trans = tmat)
    # Get transition intensities,
    Haz01 <- rbind(c(0, 0, 1),
                   ms$Haz[ms$Haz$trans == 1, ])
    Haz01$haz <- diff(c(0, Haz01$Haz))
    Haz02 <- rbind(c(0, 0, 2),
                   ms$Haz[ms$Haz$trans == 2, ])
    Haz02$haz <- diff(c(0, Haz02$Haz))
    Haz12 <- rbind(c(0, 0, 3),
                   ms$Haz[ms$Haz$trans == 3, ])
    Haz12$haz <- diff(c(0, Haz12$Haz))
    # transition probabilities from starting state at time 0
    pt <- probtrans(ms, predt = 0)[[1]] 
    
    t <- sort(unique(x$time_surv[x$event_surv == 1]))  # need to only look at events times
    t_min <- c(0, head(t, -1))
    AUC <- PT <- PI0 <- NULL
    for(i in 1:length(t)){
      # get intesity at time t[i]
      idx02 <- max(which(Haz02$time <= t[i]))
      idx12 <- max(which(Haz12$time <= t[i]))
      l02 <- Haz02$haz[idx02]
      l12 <- Haz12$haz[idx12]  
      # get transition probabilities;
      idx <- max(which(pt$time <= t[i]))
      idx_min <- max(which(pt$time <= t_min[i]))
      p01_t <- pt$pstate2[idx]
      p00_t <- pt$pstate1[idx]
      p01_tmin <- pt$pstate2[idx_min]
      p00_tmin <- pt$pstate1[idx_min]
      
      pt_est <- p01_tmin/(l02/l12*p00_tmin + p01_tmin)
      pi0_est <- p00_t/(p00_t + p01_t)
      auc_t <- 0.5 + 0.5*(pt_est - (1 - pi0_est))
      PT <- c(PT, pt_est)
      PI0 <- c(PI0, pi0_est)
      AUC <- c(AUC, auc_t)
    }
    result <- data.frame(time = t, AUC = AUC, pt = PT, pi0 = PI0)
    
    rm(x, envir = .GlobalEnv)
    return(result)
  }

# estimates incident dynamic AUC for smh model
aucID.smh <- function(data, object){
    t <- sort(unique(data$time_surv[data$event_surv == 1]))
    t_min <- c(0, head(t, -1))
    # difference between weibull and spline model
    if(is.null(dim(object$time))){
      times <- object$time
    }else{
      times <- object$time[, 1]
    }
    AUC <- PT <- PI0 <- NULL
    for(i in 1:length(t)){
      # get intesity at time t[i]
      if(length(which(times <= t[i])) != 0){
        # think about this later (hier)
        idx <- max(which(times <= t[i]))
      } else {
        idx <- min(which(times >= t[i]))
      }
      l02 <- object$intensity02[idx]
      l12 <- object$intensity12[idx]
      # get transition probabilities;
      
      # check if prediction possible
      pm <- tryCatch({predict(object, s = 0, t = t[i])$transprob}, 
                     error=function(e) e)
      if(inherits(pm, "error")) {
        i <- i - 1
        break
      }
      p01_t <- pm$Estimate[pm$Parameter == "p01"]
      p00_t <- pm$Estimate[pm$Parameter == "p00"]
      if (t_min[i] != 0){
        pm_min <- predict(object, s = 0, t = t_min[i])$transprob
        p01_tmin <- pm_min$Estimate[pm_min$Parameter == "p01"]
        p00_tmin <- pm_min$Estimate[pm_min$Parameter == "p00"]
      } else {
        p01_tmin <- 0
        p00_tmin <- 1
      }
      pt_est <- p01_tmin/(l02/l12*p00_tmin + p01_tmin)
      pi0_est <- p00_t/(p00_t + p01_t)
      auc_t <- 0.5 + 0.5*(pt_est - (1 - pi0_est))
      PT <- c(PT, pt_est)
      PI0 <- c(PI0, pi0_est)
      AUC <- c(AUC, auc_t)
    }
    result <- data.frame(time = t[1:i], AUC = AUC, pt = PT, pi0 = PI0)
    return(result)
  }


# estimates incident dynamic AUC for msm model
# data is in msm format
aucID.msm <- function(data, object, cuts, tmp_levels, marker){
    n <- length(tmp_levels)
    t <- sort(unique(data$time_surv[data$event_surv == 1]))
    t_min <- c(0, head(t, -1))
    AUC <- PT <- PI0 <- NULL
    for(i in 1:length(t)){
      # which level of time-period
      tp <- tmp_levels[which(t[i] < c(cuts, max(t) + 1) & t[i] >= c(0, cuts))]
      qm <- qmatrix.msm(object, covariates = list(1, tp), ci = "none")
      l02 <- qm["ANED", "Death"]
      l12 <- qm[marker, "Death"]
      cvs <- lapply(tmp_levels, function(x) list(1, x))
      # disease should be equal to 1
      pm <- pmatrix.piecewise.msm(object, t2 = t[i], t1 = 0,
                                  times = cpoints,
                                  covariates = cvs)
      pm_min <- pmatrix.piecewise.msm(object, t2 = t_min[i], t1 = 0,
                                      times = cpoints,
                                      covariates = cvs)
      
      p01_t <- pm["ANED", marker]
      p00_t <- pm["ANED", "ANED"]
      p01_tmin <- pm_min["ANED", marker]
      p00_tmin <- pm_min["ANED", "ANED"]
      pt_est <- p01_tmin/(l02/l12*p00_tmin + p01_tmin)
      pi0_est <- p00_t/(p00_t + p01_t)
      auc_t <- 0.5 + 0.5*(pt_est - (1 - pi0_est))
      PT <- c(PT, pt_est)
      PI0 <- c(PI0, pi0_est)
      AUC <- c(AUC, auc_t)
    }
    result <- data.frame(time = t, AUC = AUC, pt = PT, pi0 = PI0)
    return(result)
  }



