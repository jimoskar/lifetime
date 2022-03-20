### Problem 2 ----
source("~/Github/lifetime/project2/utils.R")
library(survival)
library(xtable)
library(ggplot2)
# Load dataset
tire <- read.csv("~/Github/lifetime/project2/tire.txt", sep="")

## a) ----
covariates <- c("Wedge","Inter", "Peel", "WxP")
cox.reg <- coxph(Surv(Survival, Censor) ~ Wedge + Inter + Peel + WxP, data = tire)
s <- summary(cox.reg)
A.hat <- breslow(tire, covariates, cox.reg$coefficients)

# Int. hazard rate of Weibull distr.
A.Weib <- function(t, b, k){
  return(b/k * t^k)
}

# Plot
A.hat.fn <- stepfun(A.hat$t, c(0, A.hat$a))
plot(A.hat.fn)
lines(A.hat$t, A.Weib(A.hat$t, 4e10, 15))

b <- 4e10
k <- 15
b.vec <- b * exp(as.matrix(tire[,covariates]) %*% cox.reg$coefficients)
lambda.vec <- (k/b.vec)^(1/k)
sim.T <- rweibull(34, k, lambda.vec)



mean <- lambda.vec*gamma(1+1/k) # mean of Weibull
sum(mean > 1.2)


sum(tire$Censor == 0)/max(tire$Survival)
sum(tire$Censor == 0)

sim.S <- rexp(34, rate = 1/1.2)
mean(sim.S)
sum(sim.S < sim.T)


lambda.exp <- 1.2

simulate.df <- function(df, model){
  n <- nrow(df)
  covariates <- names(model$coefficients)
  
  # Simulate survival times
  b <- 4e10
  k <- 15
  b.vec <- b * exp(as.matrix(df[,covariates]) %*% model$coefficients)
  lambda.vec <- (k/b.vec)^(1/k)
  sim.T <- rweibull(n, k, lambda.vec)
  
  # Simulate censoring times
  lambda.exp <- 1/1.2
  sim.C <- rexp(n, rate = lambda.exp)
  
  # Create simulated dataset
  censor <- (sim.C < sim.T)
  sim.df <- data.frame(df)
  sim.df$Survival <- sim.T
  sim.df$Survival[censor] <- sim.C[censor]
  sim.df$Censor <-  as.numeric(censor == 0)
  
  return(list(sim.df = sim.df, b.vec = b.vec, k = k))
}

## b) ----

# Function which performs the elimination procedure
elim <- function(df){
  y <- "Surv(Survival, Censor)"
  covariates <- c("Age", "Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")
  p.above <- TRUE
  while(p.above){
    formula <- as.formula(paste(y, paste(covariates, collapse=" + "), 
                                sep=" ~ ")) # Create formula
    model <- coxph(formula, data = df)
    s <- summary(model)
    p.vals <- s$coefficients[,5]
    if(sum(is.na(p.vals)) > 0){ # Remove covariates with NaN-values
      remove <- which(is.na(p.vals))
      covariates <- covariates[-remove]
    } else if(sum(p.vals > 0.05) == 0 || length(covariates) == 1){ # Stop elim
      p.above = FALSE
    } else{
      remove <- which(p.vals == max(p.vals))
      covariates <- covariates[-remove]
    }
  }
  return(model)
}

sim <-  simulate.df(tire, cox.reg2)
sim.df <- sim$sim.df
mod.elim <- elim(sim.df)

# Compare estimated relative risk function with the true rel. risk functions
r.est <- exp(as.matrix(tire[, names(mod.elim$coefficients)]) %*%
            mod.elim$coefficients)
r.true <- exp(as.matrix(tire[, names(cox.reg$coefficients)]) %*% cox.reg$coefficients)

log.r.df <- data.frame(log.r.est = log(r.est), log.r.true = log(r.true), x = 1:34)
ggplot(log.r.df, aes(x = x)) + geom_point(aes(y = log.r.est, color = "Estimated r")) + 
  geom_point(aes(y = log.r.true, color = "True r")) + xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("Estimated r" = "blue", "True r" = "red")) + 
  theme_minimal()

# Compare estimated integrated hazard rate with true integrated hazard rate
A0.hat <- breslow(sim.df, names(mod.elim$coefficients), mod.elim$coefficients)
r.vec <- exp(as.matrix(sim.df[, names(mod.elim$coefficients)]) %*% mod.elim$coefficients)
A.hat.df <- data.frame(a.hat1 = c(0, A0.hat$a * r.vec[1]), a.hat2 = c(0, A0.hat$a * r.vec[2]),
                       a.hat3 = c(0, A0.hat$a * r.vec[3]), t = c(0, A0.hat$t))

A.true <- function(x, idx, b.vec, k){
  return(b.vec[idx]/k * x^k)
}

# Plot for one of the cases

ggplot(A.hat.df, aes(x = t)) + geom_step(aes(y = a.hat2, color = "A.hat")) +
 geom_function(fun = A.true, args = list(idx = 2, b.vec = sim$b.vec, k = sim$k), aes(color ="A")) +
  scale_color_manual(name = " ", values = c("A.hat" = "blue", "A" = "red"), 
                     labels = expression(hat(A)(t), A(t))) + 
  coord_cartesian(ylim = c(0, 1.2)) + xlab("Study time, t") + ylab(" ") + theme_minimal()
  
## c) ----
repeat.sim <- function(num.iter, plot=FALSE){
  all.covariates <- c("Age", "Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")
  counts <- rep(0, 7)
  p <- ggplot(environment = environment())
  for(i in 1:num.iter){
    sim.df <- simulate.df(tire, cox.reg)$sim.df
    mod <- elim(sim.df)
    mod.covariates <- names(mod$coefficients)
    included <- all.covariates %in% mod.covariates
    counts[included] <- counts[included] + 1
    if(plot){
      b <- breslow(sim.df, mod.covariates, mod$coefficients)
      p <- p +  geom_step(data = data.frame(x = b$t, y = b$a),
                          aes(x = x, y = y), color = "blue", alpha = 0.5)
    }
  }
  return(list(p = p, counts = counts))
}

# Run simulations
num.iter <- 1000
c <- repeat.sim(num.iter, plot = TRUE)

# Plot proportions
ggplot(data= data.frame(prop = c$counts/num.iter, 
                        x = c("Age", "Wedge","Inter", "EB2B", "Peel", "Carbon", "WxP")), 
       aes(x=x, y=prop)) + geom_bar(stat="identity") + xlab(" ") + ylab("Proportion") + 
  theme_minimal()

# Plot baseline hazard rates
A0 <- function(x, b, k){
  return(b/k * x^k)
}
c$p + geom_function(fun = A0, args = list(b = 4e10, k = 15), color = "red", size = 1) + 
  scale_y_continuous(trans = "log") + ylab("log of integrated baseline hazard rates") +
  xlab("Study time, t") + theme_minimal()
