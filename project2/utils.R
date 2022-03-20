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