### Problem 1 ----
library(survival)
data(pbc,package="survival")
names(pbc)

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
      p.vals <- s$table[2:length(covariates), 4]
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

l <- elim(train)
l$removed
mod.elim <- l$model

r.est <- exp(as.matrix(train[, names(mod.elim$coefficients)]) %*%
               mod.elim$coefficients)

r.df <- data.frame(log.r.est = log(r.est), r.est = r.est, x = 1:312)
ggplot(r.df, aes(x = x)) + geom_point(aes(y = log.r.est, color = "Estimated r")) +
  xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("Estimated r" = "blue")) + 
  theme_minimal()


mod0 <- coxph(Surv(time, censor) ~ age + sex + ascites + hepato + spiders + edema + bili +
                          albumin + protime + alk.phos + ast + stage + trt,
      data = pbc)
sum0 <- summary(mod)
sum0$coefficients
max(sum0$coefficients[, 5])



## b) ----
lw <- elim(train, model.type = "surv")
surv.mod <- lw$model
gamma <- surv.mod$coefficients
beta <- -gamma[2:length(gamma)]/surv.mod$scale
beta
mod.elim$coefficients
r.surv <- exp(as.matrix(train[, names(beta)]) %*%
               beta)

# Plot errf together
r.df$r.surv <- r.surv
r.df$log.r.surv <- log(r.surv)
ggplot(r.df, aes(x = x)) + geom_point(aes(y = log.r.est, color = "Estimated r"))+
  geom_point(aes(y = log.r.surv, color = "Estimated r (Weibull)")) + 
  xlab("i") + ylab("log(r)") + 
  scale_color_manual(name = " ", values = c("Estimated r" = "blue", 
                                            "Estimated r (Weibull)" = "red")) + 
  theme_minimal()

