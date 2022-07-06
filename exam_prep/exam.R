# Exam TMA4275, V22
df <- data.frame(time = c(1.11, 1.13, 1.35, 1.40, 1.83, 2.16, 2.22, 2.25, 2.40, 2.65),
                 status = c(1,0,1,0,0,1,1,0,1,0))

Nelson.Aalen <- function(df, alpha, conf.int.type){
  
  # df: data frame with survival times
  # alpha: significance level of conf.int
  # conf.int.type: type of confidence interval, (1) regular or (2) log-type.
  # t: time grid
  #browser()
  times <- df$time
  n <- length(times) # Number of individuals
  status <- df$status # True if not censored
  m <- sum(status) # Number of uncensored times
  event.times <- sort(times[status == 1]) # Event times
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

ret <- Nelson.Aalen(df, 0.05, 2)
A.hat <- ret$A.hat
plot(A.hat, ylim = c(0, 3.5), xlim = c(0, 3))
plot(ret$conf.int.lower, add = TRUE, lty = "dashed")
plot(ret$conf.int.upper, add = TRUE, lty = "dashed")
  
  
