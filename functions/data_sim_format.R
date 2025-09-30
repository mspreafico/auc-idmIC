library(msm)
library(mstate)
library(SmoothHazard)

##############################
# function to simulate data
##############################

weib.simu <- function(n_sets, n_pt, adm.cens, ft,
                      a01 = 0.05, a02 = 0.05, a12 = 0.56, k = 0.5){
  
  # n_sets number of data sets of size n_pt
  # n_pt is number of patients in data
  # cens is the type of censoring
  # ft
  # a01, a02,a12 is scale of transition 01,02,12 respectively
  # k is common shape parameter of transition hazard 
  
  # generates data without censoring
  sim.data <- simu(n_sets, n_pt, a01, a02, a12, k)
  
  # censor data according to follow-up scheme
  y <- 120
  times <- seq(0, y, by = ft)
  # censore death time
  if(adm.cens){
    sim.data$event_surv <- as.numeric(sim.data$State_3 <= y)
    sim.data$time_surv <- pmin(sim.data$State_3, y)
  } else {
    ctime <- runif(nrow(sim.data), min = 5*12, max = 10*12)
    sim.data$event_surv <- as.numeric(sim.data$State_3 <= ctime)
    sim.data$time_surv <- pmin(sim.data$State_3, ctime)
  }
  # censor disease state
  obs_time <- NA
  obs_time[!is.na(sim.data$State_2)] <- sapply(sim.data$State_2[!is.na(sim.data$State_2)], 
                                               function(x) ifelse(x <= 120, 
                                                                  times[min(which(x <= times))],
                                                                  NA))
  sim.data$event_dis <- as.numeric(obs_time <= sim.data$time_surv &
                                     !is.na(obs_time))
  sim.data$time_dis <- NA
  sim.data$time_dis[sim.data$event_dis == 1] <- obs_time[sim.data$event_dis == 1]
  
  return(sim.data)
}


# samples survival times from Weibull illness-death-model without censoring
simu <- function(n_sets, n_pt, a01, a02, a12, k){
  
  # n_sets number of data sets of size n_pt
  # n_pt is number of patients in data
  # a01, a02,a12 is scale of transition 01,02,12 respectively
  # k is common shape parameter of transition hazard 
  
  n <- n_sets*n_pt
  State_1 <- numeric(n)
  
  # transition from state 1
  u <- runif(n)
  t1 <- (-log(u)/(a01 + a02))^(1/k)
  s1 <- sample(c(2, 3), n, p = c(a01/(a01 + a02), a02/(a01 + a02)), replace = TRUE)
  
  State_2 <- rep(NA, n)
  State_2[s1 == 2] <- t1[s1 == 2]
  State_2 <- State_2
  State_3 <- rep(NA, n)
  State_3[s1 == 3] <- t1[s1 == 3]
  
  # transition from state 2
  u <- runif(n)
  State_3 <- ifelse(s1 ==2 , (t1^k - log(u)/a12)^(1/k), State_3)
  State_3 <- State_3
  
  data <- data.frame(State_1 = State_1, State_2 = State_2, State_3 = State_3)
  data <- data.frame(id = 1:nrow(data), data)
  
  return(data)
  
}

##############################
# functions to format data
##############################

# long format for time-dependent Cox model
coxtd.format <- function(sim.data, n_sets, n_pt){
  
  longData <- tmerge(sim.data, sim.data, id = id,
                     death = event(time_surv, event_surv),
                     LR = tdc(time_dis))
  longData <- longData[longData$tstart < longData$tstop, ]
  # split longData into different data sets
  ids <- unique(longData$id)
  id_group <- factor(rep(c(1:n_sets), each = n_pt))
  longData$id_group <- id_group[match(longData$id, ids)]
  coxtd.data <- split(longData, longData$id_group)
  
  return(coxtd.data)
}





# long format for mstate model
mstate.format <- function(sim.data, n_sets, n_pt){
  
  idx <- is.na(sim.data$time_dis)  # fill NA values in time_dis
  sim.data$time_dis[idx] <- sim.data$time_surv[idx]
  tmat <- mstate::transMat(list(c(2, 3), c(3), c()), names = c("ANED", "LR", "Death"))
  mstate.data <- mstate::msprep(time = c(NA, "time_dis", "time_surv"),
                                status = c(NA, "event_dis", "event_surv"),
                                data = sim.data,
                                trans = tmat,
                                keep =  c("State_1", "State_2", "State_3", "time_dis", 
                                          "event_dis", "time_surv", "event_surv"))
  mstate.data <- mstate.data[mstate.data$Tstart < mstate.data$Tstop,]
  # add variabels for proportional baseline hazards
  mstate.data$disease <- ifelse(mstate.data$trans == 3, 1, 0)
  mstate.data$strata <- ifelse(mstate.data$to == 3, 2, 1)
  # split mstate.data into different data sets
  ids <- unique(mstate.data$id)
  id_group <- factor(rep(c(1:n_sets), each = n_pt))
  mstate.data$id_group <- id_group[match(mstate.data$id, ids)]
  mstate.data <- split(mstate.data, mstate.data$id_group)
  
  return(mstate.data)
}


# long format for msm
msm.format <- function(sim.data, n_sets, n_pt, ft){
  
  times <- seq(0, 120, by = ft)
  
  n <- nrow(sim.data)
  rws <- sapply(sim.data$time_surv, function(x) 1 + sum(x > times))
  id <- rep(1:n, rws)
  msm.data <- sim.data[id, ]
  msm.data$Time <- unlist(sapply(1:n, function(i) c(times[1:(rws[i] - 1)], sim.data$time_surv[i])))
  msm.data$State <- NA
  state1 <- (msm.data$event_dis == 0) |
    (msm.data$Time < msm.data$time_dis)
  state2 <- (msm.data$event_dis == 1) &
    (msm.data$Time >= msm.data$time_dis)
  state3 <- (msm.data$event_surv == 1) &
    (msm.data$Time == msm.data$time_surv)
  msm.data$State[state1] <- 1
  msm.data$State[state2] <- 2
  msm.data$State[state3] <- 3
  msm.data$disease <- ifelse(msm.data$State == 2 |
                               (msm.data$State == 3 & msm.data$event_dis == 1), 1, 0)
  # split time into two time periods
  cutpoints <- c(6, 30, 60, 90)
  msm.data$timeperiod <- ifelse(msm.data$Time < cutpoints[1], 1,
                                ifelse(msm.data$Time < cutpoints[2], 2, 
                                       ifelse(msm.data$Time < cutpoints[3], 3,
                                              ifelse(msm.data$Time < cutpoints[4], 4, 5))))
  lbls <- c(paste("[-Inf, ", cutpoints[1], ")", sep=""),
            paste("[", cutpoints[1], ", ", cutpoints[2], ")", sep = ""),
            paste("[", cutpoints[2], ", ", cutpoints[3], ")", sep = ""),
            paste("[", cutpoints[3], ", ", cutpoints[4], ")", sep = ""),
            paste("[", cutpoints[4], ", Inf)", sep = ""))
  lbls <- lbls[table(msm.data$timeperiod) > 0]
  msm.data$timeperiod <- factor(msm.data$timeperiod,
                                labels = lbls)
  
  # create list of data sets
  ids <- unique(msm.data$id)
  id_group <- factor(rep(c(1:n_sets), each = n_pt))
  msm.data$id_group <- id_group[match(msm.data$id, ids)]
  msm.data <- split(msm.data, msm.data$id_group)
  
  return(msm.data)
}


# format for SmoothHazard package
smh.format <- function(sim.data, n_sets, n_pt, ft){
  
  times <- seq(0, 120, by = ft)
  
  # create time variables  
  tt <- sim.data$time_surv
  R <- sim.data$time_dis
  R[is.na(R)] <- sapply(sim.data$time_surv[is.na(R)], 
                        function(x) times[max(which(times <= x))])
  L <- R
  L[sim.data$event_dis == 1] <- L[sim.data$event_dis == 1] - ft
  # make data set
  smh.data <- sim.data
  smh.data$d1 <- sim.data$event_dis
  smh.data$d2 <- sim.data$event_surv
  smh.data$L <- L
  smh.data$R <- R
  smh.data$tt <- tt
  # create list of data sets
  ids <- unique(smh.data$id)
  smh.data$id_group <- factor(rep(c(1:n_sets), each = n_pt))
  smh.data <- split(smh.data, smh.data$id_group)
  
  return(smh.data)
}

