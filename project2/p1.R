### Problem 1 ----
library(survival)
library(xtable)
library(ggplot2)
## a) ----
# Load dataset
tire <- read.csv("~/Github/lifetime/project2/tire.txt", sep="")

# Fit the model 1
cox.reg1 <- coxph(Surv(Survival, Censor) ~ ., data = tire)

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

covariates1 <- c("Age","Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")
covariates2 <- c("Wedge","Inter", "Peel", "WxP")
r1 <- exp(as.matrix(tire[, covariates1]) %*%
            beta.hat1)
r2 <- exp(as.matrix(tire[, covariates2]) %*% beta.hat2)

log.r.df <- data.frame(log.r1 = log(r1), log.r2 = log(r2), x = 1:34)
ggplot(log.r.df, aes(x = x)) + geom_point(aes(y = log.r1, color = "model 1")) + 
  geom_point(aes(y = log.r2, color = "model 2")) + xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("model 1" = "blue", "model 2" = "red")) + theme_minimal()

## b) ----

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

# Calculate Breslow estimator for model 1 and model 2 and plot them
b1 <- breslow(tire, covariates1, beta.hat1)
b2 <- breslow(tire, covariates2, beta.hat2)
ggplot(b1) + geom_step(aes(x = t, y = a, color = "model 1")) + 
  geom_step(data = b2, aes(x = t, y = a, color = "model 2")) +
  scale_color_manual(name = " ", values = c("model 1" = "blue", "model 2" = "red")) + 
  scale_y_continuous(trans = "log") + xlab("Study time, t") + 
  ylab(expression(log(hat(A)[0](t)))) + theme_minimal()

plot(b1, ylim = c(0,max(environment(b)$y)))

## c) Model selection ----

# Strategy: Remove covariate with lowest p-value. Stop when all covariates have p < 0.05.
summary(cox.reg1) # Remove Age

mod1 <- coxph(Surv(Survival, Censor) ~ Wedge + Inter + EB2B + Peel + Carbon + WxP, data = tire)
s <- summary(mod1) # Remove E2B2
x <- s$coefficients[, c(1,3,4,5)]
xtable(x, digits = 4)


mod2 <- coxph(Surv(Survival, Censor) ~ Wedge + Inter +  Peel + Carbon + WxP, data = tire)
s <- summary(mod2) # Remove Carbon
x <- s$coefficients[, c(1,3,4,5)]
xtable(x, digits = 4)

mod3 <- coxph(Surv(Survival, Censor) ~ Wedge + Inter +  Peel + WxP, data = tire)
summary(mod3) # All p < 0.05


covariates.mod2 <- c("Wedge","Inter", "Peel", "Carbon", "WxP")
beta.mod2 <- mod2$coefficients
r.mod2 <- exp(as.matrix(tire[, covariates.mod2]) %*%
            beta.mod2)
r2 <- exp(as.matrix(tire[, covariates2]) %*% beta.hat2)

log.r.df <- data.frame(log.r1 = log(r1), log.r.mod2 = log(r.mod2), x = 1:34)
ggplot(log.r.df, aes(x = x)) + geom_point(aes(y = log.r1, color = "model 1")) + 
  geom_point(aes(y = log.r.mod2, color = "model 2")) + xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("model 1" = "blue", "model 2" = "red")) + theme_minimal()

