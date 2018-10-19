library(raster)
library(rgdal)
library(rgeos)
library(Rcpp)
library(TMB)  
library(mvtnorm)

# eco <- readOGR("Z:/Madagascar-NMCP/Madagascar Shapefiles/Madagascar_Admin_shapefiles/Ecozones8.shp")
# 
# md <- readOGR("Z:/Madagascar-NMCP/Madagascar Shapefiles/Madagascar_Admin_shapefiles/Malareo_District.shp")
# 
# 
# W <- gTouches(md, byid = TRUE) * 1
# 
# W_eco <- gTouches(eco, byid = TRUE) * 1

setwd("C:/Users/scro3122/Documents/ModelTest/")
compile("car.cpp")
dyn.load(dynlib("car"))


N.rep <- 10
par.true <- c()
par.pred <- c()

#W_use <- W_eco
W_use <- W

##get eigenvalues
ev <- eigen(W_use)$values
l1 <- max(ev)
ln <- min(ev)
stopifnot(l1 * ln < 0)

for(i in 1:N.rep){
  print(i)
  
  N.obs <- dim(W_use)[1]
  N.beta <- 2
  phi <- runif(1, 0.1, 1/l1)
  #sigma <- runif(1, 0.1, 0.6)
  sigma <- runif(1)
  beta <- runif(N.beta)
  X <- matrix(runif(N.obs * N.beta), nrow = N.obs)
  C <- W_use * phi
  S <- solve(diag(N.obs) - C) %*% (sigma * diag(N.obs));
  

  Y <- as.vector(rmvnorm(1, mean = as.vector(X %*% beta), sigma = S))
  

  m <- MakeADFun(
    data = list(Y=runif(N.obs), X=X, Y=Y, W=W_use, phi_min=1/ln, phi_max=1/l1),
    parameters = list(sigma=runif(1), phi=runif(1), beta=runif(N.beta)),
    DLL = "car"
  )
  
  
  fit.new <- nlminb(m$par,m$fn,m$gr)
  
  # print(fit.new)
  # print(fit.new)
  # print(fit.new$par)
  # print(c(sigma, phi, beta))
  
  par.true <- rbind(par.true, c(sigma, phi, beta))
  par.pred <- rbind(par.pred, fit.new$par)
}


for(k in 1:length(fit.new$par)){
  plot(par.pred[, k], par.true[, k])
  print(cor(par.pred[, k], par.true[, k]))
}

