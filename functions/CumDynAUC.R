library(mstate)
library(msm)
library(SmoothHazard)

# Estimates cumulative dynamic AUC for mstate model
aucCD.mstate <- function(data, object, tmat, s, w){
  assign("x", data, envir = .GlobalEnv)
  ms <- msfit(object, newdata = data.frame(disease = c(0, 0, 1),
                                           trans = 1:3,
                                           strata = c(1, 2, 2)),
              trans = tmat)
  t <- s + w
  # get transition probabilities from time 0
  pm_t <- probtrans(ms, predt = 0)
  AUC <- PST <- PI0 <- NULL
  for(i in 1:length(s)){
    # get transition probabilities from time s[i]
    pm_st <- tryCatch({probtrans(ms, predt = s[i]) }, 
                      error=function(e) e)
    if(inherits(pm_st, "error")) {
      i <- i - 1
      break
    }
    # indice indicating time t and s
    idx_st <- max(which(pm_st[[1]]$time <= t[i]))
    idx_t <- max(which(pm_t[[1]]$time <= t[i]))
    idx_s <- max(which(pm_t[[1]]$time <= s[i]))
    
    p12_st <- pm_st[[2]]$pstate3[idx_st]#check
    p02_st <- pm_st[[1]]$pstate3[idx_st]#check
    p11_st <- pm_st[[2]]$pstate2[idx_st]#check
    p01_t <- pm_t[[1]]$pstate2[idx_t]
    p00_t <- pm_t[[1]]$pstate1[idx_t]
    p01_s <- pm_t[[1]]$pstate2[idx_s]
    p00_s <- pm_t[[1]]$pstate1[idx_s]
    # stop if a probability above 1
    if(sum(c(p12_st, p02_st, p11_st, p01_t,
             p00_t, p01_s, p00_s) > 1) > 0){
      i <- i - 1
      break
    }
    
    pi1_est <- p01_s*p11_st/(p00_t + p01_t)
    pi0_est <- 1 - pi1_est
    pst_est <- p12_st*p01_s/(p02_st*p00_s + p12_st*p01_s)
    
    auc_st <- 0.5 + 0.5*(pst_est - pi1_est)
    PST <- c(PST, pst_est)
    PI0 <- c(PI0, pi0_est)
    AUC <- c(AUC, auc_st)
  }
  result <- data.frame(time = s[1:i], AUC = AUC, pst = PST, pi0 = PI0)
  
  rm(x, envir = .GlobalEnv)
  return(result)
}


# estimates cumulative dynamic AUC for smh model
aucCD.smh <- function(object, s, w){
  t <- s + w
  AUC <- PST <- PI0 <- NULL
  for(i in 1:length(s)){
    # compute transition probabilities
    # check if prediction possible
    pm_st <- tryCatch({predict(object, s = s[i], t = t[i])$transprob}, 
                      error=function(e) e)
    if(inherits(pm_st, "error")) {
      i <- i - 1
      break
    }
    pm_t <- predict(object, s = 0, t = t[i])$transprob
    
    # extract transition probabilities
    p12_st <- pm_st$Estimate[pm_st$Parameter == "p12"]
    p02_st <- pm_st$Estimate[pm_st$Parameter == "p02"]
    p11_st <- pm_st$Estimate[pm_st$Parameter == "p11"]
    p01_t <- pm_t$Estimate[pm_t$Parameter == "p01"]
    p00_t <- pm_t$Estimate[pm_t$Parameter == "p00"]
    
    if (s[i] != 0){
      pm_s <- predict(object, s = 0, t = s[i])$transprob
      p01_s <- pm_s$Estimate[pm_s$Parameter == "p01"]
      p00_s <- pm_s$Estimate[pm_s$Parameter == "p00"]
    } else {
      p00_s <- 1
      p01_s <- 0
    }
    
    pi1_est <- p01_s*p11_st/(p00_t + p01_t)
    pi0_est <- 1 - pi1_est
    pst_est <- p12_st*p01_s/(p02_st*p00_s + p12_st*p01_s)
    
    auc_st <- 0.5 + 0.5*(pst_est - pi1_est)
    PST <- c(PST, pst_est)
    PI0 <- c(PI0, pi0_est)
    AUC <- c(AUC, auc_st)
  }
  if(i == 0){
    result <- NA
  } else {
    result <- data.frame(time = s[1:i], AUC = AUC, pst = PST, pi0 = PI0)
  }
  return(result)
}


# estimates cumulative dynamic AUC for msm model
# data could be in wide or long format
aucCD.msm <- function(object, ss, ww, tmp_levels, marker){
  cvs <- lapply(tmp_levels, function(x) list(1, x))
  t <- ss + ww
  AUC <- PST<- PI0 <- NULL
  for(i in 1:length(ss)){
    pm_st <- pmatrix.piecewise.msm(object, t1 = ss[i], t2 = t[i], 
                                   times = cpoints,
                                   covariates = cvs)
    pm_t <- pmatrix.piecewise.msm(object, t1 = 0, t2 = t[i], 
                                  times = cpoints,
                                  covariates = cvs)
    pm_s <- pmatrix.piecewise.msm(object, t1 = 0, t2 = ss[i], 
                                  times = cpoints,
                                  covariates = cvs)
    
    p12_st <- pm_st[marker, "Death"]
    p02_st <- pm_st["ANED", "Death"]
    p11_st <- pm_st[marker, marker]
    
    p01_t <- pm_t["ANED", marker]
    p00_t <- pm_t["ANED", "ANED"]
    
    p01_s <- pm_s["ANED", marker]
    p00_s <- pm_s["ANED", "ANED"]
    
    pi1_est <- p01_s*p11_st/(p00_t + p01_t)
    pi0_est <- 1 - pi1_est
    pst_est <- p12_st*p01_s/(p02_st*p00_s + p12_st*p01_s)
    
    auc_st <- 0.5 + 0.5*(pst_est - pi1_est)
    PST <- c(PST, pst_est)
    PI0 <- c(PI0, pi0_est)
    AUC <- c(AUC, auc_st)
  }
  result <- data.frame(time = ss, AUC = AUC, pst = PST, pi0 = PI0)
  return(result)
}

