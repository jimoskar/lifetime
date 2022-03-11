### Problem 1 ----
library(survival)
library(xtable)
library(ggplot2)
## a) ----
# Load dataset
tire <- read.csv("~/Github/lifetime/project2/tire.txt", sep="")

# Fit the model 1
covariates1 <- c("Age","Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")
cox.reg1 <- coxph(Surv(Survival, Censor) ~ ., data = tire)

cox.reg1$coefficients

# Create table for latex
s <- summary(cox.reg1)
x <- s$coefficients[, c(1,3,4,5)]
xtable(x, digits = 4)


# Fit model 2
cox.reg2 <- coxph(Surv(Survival, Censor) ~ Wedge + Inter + Peel + WxP, data = tire)

# Create table for coefficients of model 2
sum2 <- summary(cox.reg2)
sum2$coefficients
xtable(sum2$coefficients[, c(1,3,4,5)], digits = 4)


# plot r(xi, beta) for both models on log scale
beta.hat1 <- cox.reg1$coefficients
beta.hat2 <- cox.reg2$coefficients

r1 <- exp(as.matrix(tire[, c("Age","Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")]) %*%
            beta.hat1)
r2 <- exp(as.matrix(tire[, c(2, 3, 5, 7)]) %*% beta.hat2)

log.r.df <- data.frame(log.r1 = log(r1), log.r2 = log(r2), x = 1:34)
ggplot(log.r.df, aes(x = x)) + geom_point(aes(y = log.r1, color = "model 1")) + 
  geom_point(aes(y = log.r2, color = "model 2")) + xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("model 1" = "blue", "model 2" = "red")) + theme_minimal()

## b) ----

# The Breslow estimator
breslow <- function(df, covariates, beta.vec){
  x.mat <- df[, covariates]
  T.vec <- sort(df[df$Censor == 1, ]$Survival) # event times
  A.vec <- rep(NA, length(T.vec))
  A.vec[0] <- 
  for(i in 1:length(T.vec)){
    A.val <- NA
    if(i == 1){
      A.val <- 0
    } else{
      A.val <- A.vec[i - 1]
    }
    t <- T.vec[i]
    denom <- 0
    for(j in nrow(x.mat)){
      if(tire$Survival <= t){
        denom <-  denom + exp(t(beta.vec) %*% unlist(x.mat[j, ]))
      }
    }
    A.val <- A.val + 1/denom
    A.vec[i] <- A.val
  }
  return(stepfun(T.vec, A.val))
}


beta1 <- cox.reg1$coefficients
breslow(tire, covariates1, beta1)
options(error = recover)

