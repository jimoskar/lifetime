# The Breslow estimator
breslow <- function(df, covariates, beta.vec){
  x.mat <- df[, covariates] # Matrix for covariates
  r.vec <- exp(as.matrix(x.mat) %*% beta.vec) # Vector of rel. risk function estimates
  max.r <- max(r.vec)
  r.vec.norm <- r.vec/max.r # Normalize the rel. risk functions to avoid numerical problems
  T.vec <- sort(df[df$Censor == 1, ]$Survival) # Event times
  A.vec <- rep(NA, length(T.vec)) # Estimated integrated baseline hazard
  
  for(i in 1:length(T.vec)){
    A.val <- NA
    if(i == 1){
      A.val <- 0
    } else{
      A.val <- A.vec[i - 1]
    }
    t <- T.vec[i]
    idx <- tire$Survival <= t
    denom <- sum(r.vec.norm[idx])
    A.val <- A.val + 1/denom
    A.vec[i] <- A.val
  }
  
  A.vec <- A.vec/max.r # Fix scale
  return(data.frame(t = T.vec, a = A.vec))
}

Nelson.Aalen <- function(df, alpha, conf.int.type){
  
  # df: data frame with survival times
  # alpha: significance level of conf.int
  # conf.int.type: type of confidence interval, (1) regular or (2) log-type.
  # t: time grid
  
  times <- df$T.tilde
  n <- length(times) # Number of individuals
  status <- !(df$D) # True if not censored
  m <- sum(status) # Number of uncensored times
  event.times <- sort(times[status]) # Event times
  y <- rep(n, m)
  for(i in 1:m){
    y[i] <- sum(times >= event.times[i])
  }
  A.hat <- rep(0, m + 1)
  A.hat[2:(m + 1)] <- cumsum(1/y)
  sigma.hat <- rep(0, m + 1)
  sigma.hat[2:(m + 1)] = sqrt(cumsum(1/y^2))
  z = qnorm(1 - alpha/2)
  upper <- lower <- NA
  if(conf.int.type == 1){
    upper <-  A.hat + z*sigma.hat
    lower <- A.hat - z*sigma.hat
  }
  else if(conf.int.type == 2){
    exponent <- z*sigma.hat/A.hat
    exponent[is.na(exponent)] <- 0
    upper <- A.hat * exp(exponent)
    lower <- A.hat * exp(-exponent)
  }
  else{
    stop("Invalid conf.int.type")
  }
  if(length(event.times) < 1){ # No event times (for simulation)
    event.times <- 0
    A.hat <- lower <- upper <- c(0,0)}
  return(list(A.hat = stepfun(event.times, A.hat),
              conf.int.lower = stepfun(event.times, lower),
              conf.int.upper = stepfun(event.times, upper)))
}