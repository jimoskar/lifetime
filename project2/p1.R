### Problem 1 ----
## a) ----
# Load dataset
tire <- read.csv("~/Github/lifetime/project2/tire.txt", sep="")

# Fit the model 1
cox.reg1 <- coxph(Surv(Survival, Censor) ~ ., data = tire)

# Create table for latex
summary(cox.reg1)
s <- summary(cox.reg1)
x <- s$coefficients[, c(1,3,4,5)]
xtable(x)


# Fit model 2
cox.reg2 <- coxph(Surv(Survival, Censor) ~ Wedge + Inter + Peel + WxP, data = tire)
sum2 <- summary(cox.reg2)
xtable(sum2$coefficients[, c(1,3,4,5)])


# plot r(xi, beta) for both models on log scale
beta.hat1 <- cox.reg1$coefficients
beta.hat2 <- cox.reg2$coefficients

r1 <- exp(as.matrix(tire[, c(1, 2, 3, 4, 5, 6, 7)]) %*% beta.hat1)
r2 <- exp(as.matrix(tire[, c(2, 3, 5, 7)]) %*% beta.hat2)

log.r.df <- data.frame(log.r1 = log(r1), log.r2 = log(r2), x = 1:34)
ggplot(log.r.df, aes(x = x)) + geom_point(aes(y = log.r1)) + 
  geom_point(aes(y = log.r2), color = "red")

