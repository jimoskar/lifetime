### Problem 1 ----
library(survival)
data(pbc,package="survival")
names(pbc)

## a) ----
# Add censor column to dataset and pick out first 312 observations
censor <- pbc$status
censor[censor == 1] = 0
censor[censor == 2] = 1
pbc$censor <- censor
pbc <- pbc[1:312, ]

mod0 <- coxph(Surv(time, censor) ~ age + sex + ascites + hepato + spiders + edema + bili +
                          albumin + protime + alk.phos + ast + stage + trt,
      data = pbc)
sum0 <- summary(mod)
sum0$coefficients
max(sum0$coefficients[, 5])

# Elimination strategy: Remove highest p-value (except for treatment) until all <= 0.05.
# 1) Remove alk.phos
mod1 <- coxph(Surv(time, censor) ~ age + sex + ascites + hepato + spiders + edema + bili +
                albumin + protime + ast + stage + trt,
              data = pbc)
sum1 <- summary(mod1)
sum1
# 2) Remove spiders
mod2 <- coxph(Surv(time, censor) ~ age + sex + ascites + hepato  + edema + bili +
                albumin + protime + ast + stage + trt,
              data = pbc)
sum2 <- summary(mod2)
sum2
# 3) Remove ascites
mod3 <- coxph(Surv(time, censor) ~ age + sex + hepato  + edema + bili +
                albumin + protime + ast + stage + trt,
              data = pbc)
sum3 <- summary(mod3)
sum3
# 4) Remove hepato
mod4 <- coxph(Surv(time, censor) ~ age + sex  + edema + bili +
                albumin + protime + ast + stage + trt,
              data = pbc)
sum4 <- summary(mod4)
sum4
# 5) Remove sex
mod5 <- coxph(Surv(time, censor) ~ age + edema + bili +
                albumin + protime + ast + stage + trt,
              data = pbc)
sum5 <- summary(mod5)
sum5

## b) ----
wmod0 <- survreg(Surv(time, censor) ~ age + sex + ascites + hepato + spiders + edema + bili +
                   albumin + protime + alk.phos + ast + stage + trt,
                 data = pbc)
summary(wmod0)
