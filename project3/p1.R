### Problem 1 ----
library(survival)
library(xtable)
library(ggplot2)
data(pbc,package="survival")

## a) ----

# Construct train dataset
train <- pbc[1:312, 
             c("age", "sex", "ascites", "hepato", "spiders", "edema",
               "bili", "albumin", "protime", "alk.phos", "ast", "stage",
               "trt", "time")] # pick out first 312 observation
# Add censor column
censor <- pbc$status[1:312]
censor[censor == 1] = 0
censor[censor == 2] = 1
train$censor <-  censor
# Add square and log() of continuous variables
train$age2 <- (train$age)^2
train$lage <- log(train$age)
train$bili2 <- (train$bili)^2
train$lbili <- log(train$bili)
train$albumin2 <- (train$albumin)^2
train$lalbumin <- log(train$albumin)
train$protime2 <- (train$protime)^2
train$lprotime <- log(train$protime)
train$alk.phos2 <- (train$alk.phos)^2
train$lalk.phos <- log(train$alk.phos)
train$ast2 <- (train$ast)^2
train$last <- log(train$ast)


# Elimination strategy: Remove highest p-value until all <= 0.05.
# Function for backward elimination of features:
elim <- function(df, model.type = "cox"){
  y <- "Surv(time, censor)"
  covariates <- setdiff(colnames(df), c("time", "censor"))
  removed <- c()
  p.above <- TRUE
  while(p.above){
    formula <- as.formula(paste(y, paste(covariates, collapse=" + "), 
                                sep=" ~ ")) # Create formula
    model <- p.vals <- NA
    if(model.type == "cox"){
      model <- coxph(formula, data = df)
      s <- summary(model)
      p.vals <- s$coefficients[,5]
    } else if(model.type == "surv"){
      model <- survreg(formula, data = df)
      s <- summary(model)
      p.vals <- s$table[2:(length(covariates)+1), 4]
    } else{
      stop("Invalid model.")
    }
    if(sum(p.vals > 0.05) == 0 || length(covariates) == 1){ # Stop elim
      p.above = FALSE
    } else{
      remove <- which(p.vals == max(p.vals))
      removed <- c(removed, covariates[remove])
      covariates <- covariates[-remove]
    }
  }
  return(list(model = model, removed = removed))
}

# Run step-wise backward elimination procedure
elim.cox <- elim(train)
elim.cox$removed
mod.cox <- elim.cox$model

r.cox <- exp(as.matrix(train[, names(mod.cox$coefficients)]) %*%
               mod.cox$coefficients)

# Plot errf
r.df <- data.frame(log.r.cox = log(r.cox), r.cox = r.cox, x = 1:312)
ggplot(r.df, aes(x = x)) + geom_point(aes(y = log.r.cox, color = "Estimated r")) +
  xlab("Case no. i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("Estimated r" = "blue")) + 
  theme_minimal()
ggsave("figures/errf_cox.pdf", height = 5, width = 8)

## b) ----
elim.surv <- elim(train, model.type = "surv")
elim.surv$removed
mod.surv <- elim.surv$model

# Find corresponding coefficients
gamma <- mod.surv$coefficients
beta <- -gamma[2:length(gamma)]/mod.surv$scale
r.surv <- exp(as.matrix(train[, names(beta)]) %*%
               beta)

# Plot errf together
r.df$r.surv <- r.surv
r.df$log.r.surv <- log(r.surv)
ggplot(r.df, aes(x = x)) + 
  geom_point(aes(y = log.r.cox, color = "Estimated r (Cox)"))+
  geom_point(aes(y = log.r.surv, color = "Estimated r (Weibull)")) + 
  xlab("Case no. i") + ylab("log(r)") + 
  scale_color_manual(name = " ", 
                     values = c("Estimated r (Cox)" = "blue", 
                                            "Estimated r (Weibull)" = "red")) + 
  theme_minimal()
ggsave("figures/errf_both.pdf", height = 5, width = 8)

## c) ----
# Construct table with coefficients
sum.surv <- summary(mod.surv)
p.vals <- sum.surv$table[2:(length(mod.surv$coefficients)),4]
df <- data.frame("coef" = beta, "p-value" = p.vals)
xtable(df, digits = 4)


## d) ----
marg.var <- diag(sum.surv$var)
z.alpha <- qnorm(0.05, lower.tail = FALSE)
params <- sum.surv$table[,1]
lower <- params - z.alpha*sqrt(marg.var)
upper <- params + z.alpha*sqrt(marg.var)
conf.ints <- data.frame(upper = upper, lower = lower)
xtable(conf.ints, digits = 5)


# Remove??
f <- function(phi, gamma, x){
  val <- exp(-phi)*t(gamma)%*%x
  return(val)
}

var.T <- function(phi, gamma, x, var.phi, var.g, cov.pg){
  g.x <- t(gamma)%*%x
  val <- exp(-2*phi)*(g.x^2*var.phi + sum(x^2*var.g) -
                        2*g.x*sum(x*cov.pg))
  return(val)
}

# Find variance of errf
var.vec <- rep(NA, 312)
covariates <- names(tail(mod.surv$coefficients,-1))
phi <- params[8]
gamma <- params[2:8]
var.phi <- marg.var[9]
var.g <- marg.var[2:8]
cov.pg <- mod.surv$var[9, 2:8]
for(i in 1:312){
  x <- train[i, covariates]
  var.vec[i] <- var.T(phi, as.numeric(gamma), as.numeric(x),
                      var.phi, var.g, cov.pg) 
}

r.df$upper <- r.df$log.r.surv + z.alpha*sqrt(var.vec)
r.df$lower <- r.df$log.r.surv - z.alpha*sqrt(var.vec)

ggplot(r.df, aes(x = x)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, colour="Confidence interval"), width=2, size = 0.5) + 
  geom_point(aes(y = log.r.surv, color = "Estimated r (Weibull)"), size = 0.2) +
  xlab("Case no. i") + ylab(" ") + 
  scale_color_manual(name = " ", 
                     values = c("Estimated r (Weibull)" = "blue4", 
                                "Confidence interval" = "slategray2" ),
                     labels = expression(hat(beta)%.%x, paste("C.I."))) + 
  theme_minimal()
ggsave("figures/errf_confint.pdf", height = 4, width = 10)



