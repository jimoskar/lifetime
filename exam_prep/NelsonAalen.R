### Solution to 3.1 in ABG

#P for placebo
#M is the 6-MP group

T_P = c(0,1,2,3,4,5,8,11,12,15,17,22,23) #Event times
d_P = c(0,2,2,1,2,2,4,2,2,1,1,1,1) #Number of events
c_P = rep(0,length(T_P)) #Censored observations

T_M = c(0,6,7,9,10,11,13,16,17,19,20,22,23,25,32,34,35)
d_M = c(0,3,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0)
c_M = c(0,1,0,1,1,1,0,0,1,1,1,0,0,1,2,1,1)

Nelson_Aalen = function(Times,d,c){
  # Calculates the Nelson-Aalen estimate according to 3.12 in ABG
  m = length(Times)
  Y = rep(21,m)
  Y[2:m] = Y[2:m] - cumsum(d+c)[1:(m-1)]
  
  DeltahatA = rep(0,m)
  for (j in 2:m){
    # This is the adjusted inner sum
    if (d[j] > 0){
      for (l in 0:(d[j]-1)){
        DeltahatA[j] = DeltahatA[j] + 1/(Y[j]-l)
      }
    }
  }
  return(stepfun(Times[2:m],cumsum(DeltahatA))) #The outer sum is here
}

hatA_P = Nelson_Aalen(T_P,d_P,c_P)
hatA_M = Nelson_Aalen(T_M,d_M,c_M)


t = seq(0,40,by=0.01)

plot(t,hatA_P(t),type="l",col="blue",ylab="Nelson-Aalen estimate")
lines(t,hatA_M(t),col="red")
legend(0,3.5,c("Placebo","6-MP"),col=c("blue","red"),lty=1)
