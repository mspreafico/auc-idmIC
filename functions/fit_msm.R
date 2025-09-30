library(msm)

##############################
# function to fit msm model
##############################

fit_msm <- function(data, 
                    method = "BFGS",
                    gen.inits = TRUE,
                    control = list(fnscale = 25000, maxit = 10000)){
  fm <- msm(State ~ Time, subject = id, data = data,
            qmatrix = qmat, death = 3, 
            covariates = ~ disease + timeperiod,
            fixedpars = c(3, 4), 
            constraint = list("timeperiod[6, 30)" = c(1, 2, 2),
                              "timeperiod[30, 60)" = c(1, 2, 2),
                              "timeperiod[60, 90)" = c(1, 2, 2),
                              "timeperiod[90, Inf)" = c(1, 2, 2)),
            qconstraint = c(1, 2, 2),
            method = method, 
            gen.inits = gen.inits,
            control = control)
  return(fm)
}
