### Solution to problem 3.7 in ABG


Times = c(0,1,2,3,4,5,8,11,12,15,17,22,23) #Event times
d = c(0,2,2,1,2,2,4,2,2,1,1,1,1) #Number of events
c = rep(0,length(Times)) #Censored observations

Kaplan_Meier = function(Times,d,c){
  # Calculates the Kaplan-Meier estimate according to 3.32 in ABG
  m = length(Times)
  Y = rep(21,m)
  Y[2:m] = Y[2:m] - cumsum(d+c)[1:(m-1)]
  
  S = cumprod(1 - d/Y)
  
  Shat = stepfun(Times[2:m],S)
  
  return(Shat)
}

SigmaHat2 = function(Times,d,c){
  # Calculates the variance estimate in 3.14/3.15,
  # which will then be used to calculate the estimated variance
  # of the Kaplan-Meier estimates
  m = length(Times)
  Y = rep(21,m)
  Y[2:m] = Y[2:m] - cumsum(d+c)[1:(m-1)]
  
  Deltahatsigma2 = rep(0,m)
  for (j in 1:m){
    # This is the adjusted inner sum
    if (d[j] > 0){
      for (l in 0:(d[j]-1)){
        Deltahatsigma2[j] = Deltahatsigma2[j] + 1/(Y[j]-l)^2
      }
    }
  }
  hatsigma2 = stepfun(Times[2:m],cumsum(Deltahatsigma2))
  
  return(hatsigma2)
}

hatS = Kaplan_Meier(Times,d,c)

t = seq(0,max(Times),by=0.01)

hatTau2 = hatS(t)^2*SigmaHat2(Times,d,c)(t)
exponent = qnorm(0.975)*sqrt(hatTau2)/(hatS(t)*log(hatS(t)))
conf_upper = hatS(t)*exp(-exponent)
conf_lower = hatS(t)*exp(exponent)

par(mfrow = c(1,1))
plot(t,hatS(t),type="l",col="blue",ylab="Kaplan-Meier estimate")
lines(t,conf_upper,lty=2,col="blue")
lines(t,conf_lower,lty=2,col="blue")

