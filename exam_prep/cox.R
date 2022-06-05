### Cox regression with pneumonia data

pneu=read.table("http://www.math.ntnu.no/emner/TMA4275/2021v/Datasets/pneumonia.txt",header=T)

## b) ---- 
# something is wrong here
pneu
sorted <- pneu[order(pneu$time),]
sorted

partlik <- function(beta){
  log.lik <- 0
  for(j in 1:nrow(pneu)){
    if(sorted$cens[j] == 0){
      next 
    }
    Tj <- sorted$time[j]
    risk.sum <- 0
    for(l in 1:nrow(pneu)){
      if(pneu$time[l] <= Tj){
        risk.sum <- risk.sum + exp(beta*pneu$x[l])
      }
    }
    log.lik <- log.lik + (beta*sorted$x[j] - log(risk.sum))
  }
  return(exp(log.lik))
}

betas = seq(-3,3,length.out = 100)
lhs = array(0,100)
for (i in 1:100){
  lhs[i] = partlik(betas[i])
}

plot(betas,lhs,type="l")
k = which.max(lhs)
lines(c(betas[1],betas[k]),c(lhs[k],lhs[k]),lty=2,col="red")
lines(c(betas[k],betas[k]),c(0,lhs[k]),lty=2,col="red")

## c) ----
library(survival)
fit.pneu = coxph(Surv(time,cens==1)~x,method="breslow",data=pneu)
s <- summary(fit.pneu)

## d) ----
# p-value is 0.0571, which means that we reject the hypothesis that there is
# a sign. difference between the groups at sign.level 0.05.

# estimated relative risk
s$coefficients[2]

## e) ----
Times = c( 2, 3, 4, 6, 9,10,11,12,17,23,24,26,32) #Times of interest (events and censored observations)
d_P =   c( 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1) #Number of events, Pneumonia
c_P =   c( 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0) #Number of censored obs, Pneu
d_M =   c( 1, 0, 0, 2, 0, 1, 1, 0, 0, 1, 0, 0, 0) #Number of events, no pneu
c_M =   c( 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0) #Number of censored obs, no pneu

d = d_P + d_M
c = c_P + c_M

#Number of individuals in each group
numP = 7 
numM = 8

m = length(Times)

#Calculate the number at risk at each time
Y_P = rep(numP,m)
Y_M = rep(numM,m)
Y_P[2:m] = Y_P[2:m] - cumsum(d_P+c_P)[1:(m-1)]
Y_M[2:m] = Y_M[2:m] - cumsum(d_M+c_M)[1:(m-1)]
Y = Y_P + Y_M

#Calculate the test statistic
Z_1 = sum(d_P*Y_M/Y - d_M*Y_P/Y)
V_11 = sum(d*Y_M*Y_P/Y^2)
X = Z_1^2/V_11

#Get the p-value
pchisq(X,1,lower.tail = F)

#Checking the expected values for each group
E_P = sum(d*Y_P/Y) #Expected number of observations in the Pneumonia group
E_M = sum(d*Y_M/Y) #Expected number of observations in the other group

## f) ----
survdiff(Surv(time,cens)~x,data=pneu)
