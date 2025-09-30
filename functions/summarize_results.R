#########################
# functions to summarize results
########################

model.perf <- function(model, true.coef, msm.type=F){
  
  coefs <- NULL
  if(msm.type){
    for(i in 1:length(model)){
      tmp <- tryCatch(model[[i]]$estimates[6], error = function(e) NA)
      coefs <- c(coefs, tmp)
    }}else{
      for(i in 1:length(model)){
        tmp <- tryCatch(model[[i]]$coefficients, error = function(e) NA)
        coefs <- c(coefs, tmp)
      }
    }
  
  result = list()
  result$valid <- sum(!is.na(coefs))
  result$LR <- mean(coefs, na.rm = TRUE)
  result$seLR <- sd(coefs, na.rm = TRUE)
  result$bias <- result$LR - true.coef
  result$mse <- (result$bias)^2 + (result$seLR)^2
  result$rmse <- sqrt(result$mse)
  result$LR_HR <- exp(result$LR)
  
  return(result)
  
}

auc.at.t <- 
  # function that computes average AUC at specific times 
  function(object.list, true.auc, time.var = "time",
           times = 0:5*12, except = NULL, inc = TRUE){
    # removes invalid estimation
    if(!is.null(except)){
      for(j in 1:length(except)){
        # sometimes transition probabilities of last value are above one
        # except contains vector of problematic objects
        if(inc){
          object.list[[except[j]]]$AUC[length(object.list[[except[j]]]$AUC)] <- NA
        } else {
          object.list[[except[j]]]$AUC[61] <- NA
        }
      }
    }
    auc <- NULL
    for(i in 1:length(object.list)){
      object <- object.list[[i]]
      if((class(object) == "try-error") | (length(object) == 1)) {
      } else {
        # computes index of time which is <= times elements
        # to find correct AUC, also checks if the time is 
        # within 6 months of x
        tt <- object[[time.var]]
        idx <- numeric(length(times))        
        for(j in 1:length(times)){
          tmp <- which(tt <= times[j])
          idx[j] <- ifelse(length(tmp) == 0 | 
                            (times[j] - tt[max(tmp)]) > 6, NA, max(tmp))
        }
        auc <- rbind(auc, ifelse(is.na(idx), NA, object$AUC[idx]))
      }
    }
    AUC <- apply(auc, 2, function(x) mean(x, na.rm = T))
    se <- apply(auc, 2, function(x) sd(x, na.rm = T))
    valid <- apply(auc, 2, function(x) sum(!is.na(x)))
    bias <- AUC - true.auc
    rmse <- sqrt(bias^2 + se^2)
    res <- data.frame(cbind(times, AUC, valid, bias, se, rmse))
    colnames(res) <- c("time", "AUC", "valid", "bias", "se", "rmse")
    return(res)
  }


plot.all <- function(object.list, time.var="time", cumulative = FALSE, 
                     title="", ylim=c(0.4,0.8), xlim = c(0, 10),
                     except = NULL, inc = TRUE){
  # removes invalid estimation
  if(!is.null(except)){
    for(j in 1:length(except)){
      # sometimes transition probabilities of last value are above one
      # except contains vector of problematic objects
      if(inc){
        object.list[[except[j]]]$AUC[length(object.list[[except[j]]]$AUC)] <- NA
      } else {
        object.list[[except[j]]]$AUC[61] <- NA
      }
    }
  }
  j <- 1
  while((class(object.list[[j]]) == "try-error") | (length(object.list[[j]]) == 1)){
    j <- j + 1
  }
  # plot lines
  plot(object.list[[j]][[time.var]]/12, object.list[[j]]$AUC, 
       xlab = "Time t in years", ylab = "AUC(t)", main = title, type = "l", 
       col = rgb(0,1,0,0.3), xlim = xlim, ylim = ylim, cex.lab = 2)
  for(i in 1:length(object.list)){
    if((class(object.list[[i]]) == "try-error") | (length(object.list[[i]]) == 1)) {
    } else {
      lines(object.list[[i]][[time.var]]/12, 
            object.list[[i]]$AUC, type = "l", 
            col = rgb(0,1,0,0.3))
    }
  }
}

find.except.mstate <- function(object.list){
  except = NULL
  for(i in 1:length(object.list)){
    nrow = dim(object.list[[i]])[1]
    if(length(which(object.list[[i]][nrow,2:4]>1))>0 | length(which(object.list[[i]][nrow,2:4]<0))>0){
      except = c(except,i)
    }
  }
  return(except)
}
