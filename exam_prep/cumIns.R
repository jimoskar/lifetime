library(survival)
library(cmprsk)

#Number of simulations
nsim <- 1000

simulate_and_plot = function(nsim,lambda1,lambda2,seed){
  #Simulate data
  set.seed(seed)
  T1s <- rexp(nsim,lambda1)
  T2s <- rexp(nsim,lambda2)
  
  obstime <- array(0,nsim)
  obstime <- pmin(T1s,T2s)
  d <- 1*(obstime==T1s)+2*(obstime==T2s)
  
  causespec1 <- survfit(Surv(obstime,d==1)~1,type="fh2")
  plot(causespec1,fun="cumhaz",mark.time=F,
       ylab=paste("Lambda1: ",lambda1,", Lambda2: ",lambda2,sep=""),ylim=c(0,lambda1*max(causespec1$time))) 
  
  #Estimate and plot the cimulative incidence functions:
  ci = cuminc(obstime,d)
  plot(ci$`1 1`$time,ci$`1 1`$est,lty=1,type="s",ylim=c(0,1),ylab="")
  lines(ci$`1 2`$time,ci$`1 2`$est,lty=2,type="s",ylab="")
  legend(0,1,c("T1","T2"),lty=1:2,bty="n")
}

par(mfrow=c(3,2))
simulate_and_plot(nsim,lambda1=2,lambda2=1,seed=4275)
simulate_and_plot(nsim,lambda1=2,lambda2=2,seed=4275)
simulate_and_plot(nsim,lambda1=2,lambda2=5,seed=4275)
