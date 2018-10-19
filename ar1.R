setwd("C:/Users/scro3122/Documents/ModelTest/")
library(TMB)


compile("ar1.cpp")
dyn.load(dynlib("ar1"))
N <- 200
phi.true <- c()
phi.pred <- c()

par.true <- c()
par.pred <- c()

set.seed(9)
ptm <- proc.time()
for(j in 1:N){
  print(j)
  x <- c()
  x[1] <- runif(1)
  phi <- runif(1, -1, 1)
  #phi.true[j] <- phi
  
  M <- 100
  sigma <- runif(1)
  for(i in 2:M){
    x[i] <- phi*x[i-1] + rnorm(1, sd = sigma)
  }
  
  
  #plot(x)
  
  
  
  m <- MakeADFun(
    data = list(x=x, sd=0.1),
    parameters = list(inv_phi=runif(1), sigma=runif(1)),
    DLL = "ar1"
  )
  m$env$tracepar <- TRUE
  
  fit.new <- nlminb(m$par,m$fn,m$gr)
  #,
  #                  control=list(iter.max=300,eval.max=300))
  #report <- sdreport(m, getJointPrecision = 1)
  
  #print(fit.new)
  
  
  #phi.pred[j] <- fit.new$par
  par.pred <- rbind(par.pred, fit.new$par)
  par.true <- rbind(par.true, c(phi, sqrt(sigma)))
  
}

for(k in 1:2){
  plot(par.true[, k], par.pred[, k])
  abline(a=0, b=1)
}
# plot(phi.true, phi.pred)
# print(cor(phi.true, phi.pred))

print(proc.time() - ptm)
