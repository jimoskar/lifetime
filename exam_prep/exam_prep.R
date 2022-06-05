### Exam problems
## 3. ----
library(survival) 
library(ggplot2)
library(ggfortify)
df <- read.table("leukemia.txt", header = TRUE)
km.placebo <- survfit(Surv(time[treat == 1], status[treat == 1]) ~ 1, data=df)
km.treat <- survfit(Surv(time, status) ~ treat, data=df)
autoplot(km.treat)

### Exercise 1.4
# Simulate Poisson process

sim.pois <- function(t, lambda){
  N <- c(0)
  jump.time <- c()
  agg.time <- 0
  stop = FALSE
  while(!stop){
    sim.exp <- rexp(1, lambda)
    N <- c(N, tail(N, 1) + 1)
    jump.time <- c(jump.time, agg.time + sim.exp)
    if(agg.time + sim.exp > t){
      stop <- TRUE
    } else{
      agg.time <- agg.time + sim.exp
    }
  }
  return(list(t = head(jump.time, -1), N = head(N, -1)))
}

l <- sim.pois(25, 1)
plot(stepfun(l$t,l$N))


