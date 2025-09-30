# Code to compute theoretical values of incident dynamic and 
# cumulative dynamic AUC and related quantities.

#########################
# hazard functions
#########################

lmd <- 
  # weibul hazard
  function(a, k, t){
    return(a*k*t^(k-1))
  }

Lmd <- 
  # cumulative weibul hazard
  function(a, k, t){
    return(a*t^k)
  }

#########################
# Transition probabilities
#########################

S0 <- 
  # probability of staying in state 0
  function(a01, a02, k, t){
    return(exp(-(a01 + a02)*t^k))
  }
S1 <- 
  # probability of staying in state 1
  function(a12, k, t){
    return(exp(-a12*t^k))
  }

P01 <- 
  # probability of being in state 1 
  # conditional on being in state 0 at time 0
  function(a01, a02, a12, k, s = 0, t){
    s1t <- S1(a12 = a12, k = k, t = t)
    s1s <- S1(a12 = a12, k = k, t = s)
    s0t <- S0(a01 = a01, a02 = a02, k = k, t = t)
    s0s <- S0(a01 = a01, a02 = a02, k = k, t = s)
    if(a01 + a02 == a12){
      p01st <- a01*(s1t/s1s*t^k - s0t/s0s*s^k)
    } else {
      p01st <- a01/(a01 + a02 - a12)*(s1t/s1s - s0t/s0s)
    }
    return(p01st)
  }

P02_0 <- 
  # probability of being in state 2 at time t 
  # conditional on being in state 0 at time s
  # and not experiencing disease (state 1)
  function(a01, a02, k, s = 0, t){
    s0t <- S0(a01 = a01, a02 = a02, k = k, t = t)
    s0s <- S0(a01 = a01, a02 = a02, k = k, t = s)
    p02_0st <- a02/(a01 + a02)*(1 - s0t/s0s)
    return(p02_0st)
  }

P02_1 <- 
  # probability of being in state 2 
  # conditional on being in state 0 at time s
  # and experienceing disease (state 1)
  function(a01, a02, a12, k, s = 0, t){
    s0t <- S0(a01 = a01, a02 = a02, k = k, t = t)
    s0s <- S0(a01 = a01, a02 = a02, k = k, t = s)
    s1t <- S1(a12 = a12, k = k, t = t)
    s1s <- S1(a12 = a12, k = k, t = s)
    part1 <- a01/(a01 + a02)*(1 - s0t/s0s)
    if(a01 + a02 == a12){
      part2 <- a01*(s1t/s1s*t^k - s0t/s0s*s^k)
    } else {
      part2 <- a01/(a01 + a02 - a12)*(s1t/s1s - s0t/s0s)
    }
    p02_1st <- part1 - part2
    return(p02_1st)
  }

P02 <- 
  # probability of being in state 2 
  # conditional on being in state 0 at time s
  function(a01, a02, a12, k, s = 0, t){
    p02_0st <- P02_0(a01 = a01, a02 = a02, k = k, s = s, t = t)
    p02_1st <- P02_1(a01 = a01, a02 = a02, a12 = a12, k = k, s = s, t = t)
    p02st <- p02_0st + p02_1st
    return(p02st)
  }

P12 <- 
  # probability of being in state 2 
  # conditional on being in state 1 at time s
  function(a12, k, s = 0, t){
    s1t <- S1(a12 = a12, k = k, t = t)
    s1s <- S1(a12 = a12, k = k, t = s)
    p12st <- 1 - s1t/s1s
    return(p12st)
  }

P00 <- 
  # probability of being in state 0 at time t 
  # conditional on being in state 0 at time s
  function(a01, a02, k, s = 0, t){
    s0t <- S0(a01 = a01, a02 = a02, k = k, t = t)
    s0s <- S0(a01 = a01, a02 = a02, k = k, t = s)
    p00st <- s0t/s0s
    return(p00st)
  }

P11 <- 
  # probability of being in state 1 at time t
  # conditional on being in state 1 at time s
  function(a12, k, s = 0, t){
    s1t <- S1(a12 = a12, k = k, t = t)
    s1s <- S1(a12 = a12, k = k, t = s)
    p11st <- s1t/s1s
    return(p11st)
  }

#########################
# functions for AUCs:
# p(t), pi_1(t), pi_0(t), p(s, t), pi_1(s, t), pi_0(s, t)
#########################
pt <- 
  # probability that a person that dies at t has history of illness
  function(a01, a02, a12, k, t) {
    lmd12t <- lmd(a = a12, k = k, t = t)
    lmd02t <- lmd(a = a02, k = k, t = t)
    p01t <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = t)
    p00t <- P00(a01 = a01, a02 = a02, k = k, t = t)
    ptt <- lmd12t*p01t/(lmd02t*p00t + lmd12t*p01t)
    return(ptt)
  }

pi1 <- 
  # prevalence of illnes at time t
  function(a01, a02, a12, k, t) {
    p01t <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = t)
    p00t <- P00(a01 = a01, a02 = a02, k = k, t = t)
    pi1t <- p01t/(p00t + p01t)
    return(pi1t)
  }

pi0 <- 
  # 1 - prevalence of illness
  function(a01, a02, a12, k, t) {
    pi1t <- pi1(a01 = a01, a02 = a02, a12 = a12, k = k, t = t)
    pi0t <- 1 - pi1t
    return(pi0t)
  }

pi0st <-  function(a01, a02, a12, k, s, t){
  p00t <- P00(a01 = a01, a02 = a02, k = k, t = t)
  p00s <- P00(a01 = a01, a02 = a02, k = k, t = s)
  p01st <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, s = s, t = t)
  p01t <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = t) 
  num <- p00t + p00s*p01st
  denom <- p00t + p01t
  return(num/denom)
}

pi1st <-  function(a01, a02, a12, k, s, t){
  p11st <- P11(a12 = a12, k = k, s = s, t = t)
  p01s <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = s)
  p00t <- P00(a01 = a01, a02 = a02, k = k, t = t)
  p01t <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = t) 
  num <- p11st*p01s
  denom <- p00t + p01t
  return(num/denom)
}

pst <-  function(a01, a02, a12, k, s, t){
  p12st <- P12(a12 = a12, k = k, s = s, t = t)
  p01s <- P01(a01 = a01, a02 = a02, a12 = a12, k = k, t = s)
  p02st <- P02(a01 = a01, a02 = a02, a12 = a12, k = k, s = s, t = t)
  p00s <- P00(a01 = a01, a02 = a02, k = k, t = s)
  num <- p12st*p01s
  denom <- p02st*p00s + p12st*p01s
  return(num/denom)
}

#########################
# AUC functions
#########################

weib.auc.id <- 
  function(a01, a02, a12, k, t){
    ptt <- pt(a01 = a01, a02 = a02, a12 = a12, k = k, t = t)
    pi1t <- pi1(a01 = a01, a02 = a02, a12 = a12, k = k, t = t)
    auct <- 0.5 + 0.5*(ptt - pi1t)
    return(auct)
  }

weib.auc.cd <- function(a01, a02, a12, k, s, t){
  p <- pst(a01 = a01, a02 = a02, a12 = a12, k = k, s = s, t = t)
  pi1 <- pi1st(a01 = a01, a02 = a02, a12 = a12, k = k, s = s, t = t)
  auct <- 0.5 + 0.5*(p - pi1)
  return(auct)
}
