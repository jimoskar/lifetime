library(survival)
set.seed(4275)

n = 100
x=rgamma(n,2)
T=sqrt(rexp(n)*2*exp(-x))
C=rexp(n,0.5)
t=pmin(C,T)
delta=1*(T<C)

cfit = coxph(Surv(t,delta)~x)
summary(cfit)
martres = cfit$residuals
plot(x,martres)
lines(lowess(x,martres),col="red",lty=2)

## c) ----
T <- sqrt(rexp(n) *2 * exp(-log(x)))
C=rexp(n,0.5)
t=pmin(C,T)
delta=1*(T<C)
cfit = coxph(Surv(t,delta)~x)
summary(cfit)
martres = cfit$residuals
plot(x,martres)
lines(lowess(x,martres),col="red",lty=2)
