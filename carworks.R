setwd("C:/Users/scro3122/Documents/ModelTest")
library(rgeos)
library(raster)
library(mvtnorm)
library(rgdal)
library(ggplot2)
library(TMB)


# W <- shapefile("C:/Users/scro3122/Documents/Madagascar_Admin_shapefiles/Malareo_Communes.shp")
# 
# admin2 <- shapefile("C:/Users/scro3122/Documents/Madagascar_Admin_shapefiles/Malareo_Region.shp")
# W <- gTouches(admin2, byid=TRUE) * 1

#get eigenvalues to bound phi
e.values <- eigen(W)$values
l1 <- max(e.values)
ln <- min(e.values)
stopifnot(l1 * ln < 0)

compile("carworks.cpp")
dyn.load(dynlib("carworks"))

# compile("car_2.cpp")
# dyn.load(dynlib("car_2"))

# compile("car_simple.cpp")
# dyn.load(dynlib("car_simple"))

par.pred <- c()
par.true <- c()

# World2 <- fortify(admin2)
# ids <- as.numeric(World2$id)
# Y.pixel <- rep(Y - X%*%beta, sapply(0:21, function(i) sum(ids == i)))
# World2$Y <- Y.pixel
# p <- ggplot() + 
#   geom_polygon(data = World2, aes(x = long, y = lat, group = group, fill=Y.pixel),
#                colour = "black", size = 0.5) + coord_fixed()
# print(p)
# print(phi)

#set.seed(10)

#seed <- sample.int(100, 1)
#set.seed(seed)
#set.seed(9)

N.reps <- 100
for(k in 1:N.reps){
  print(k)
  N.beta <- 2
  N.obs <- dim(W)[1]
  #phi <- runif(1, 1/ln, 1/l1)
  phi <- runif(1, 0.1, 1/l1)
  C <- phi * W
  variance <- runif(1)
  M <- diag(N.obs) * variance
  X <- matrix(runif(N.beta * N.obs), ncol=N.beta)
  beta <- runif(N.beta) 
  Y <- as.vector(rmvnorm(1, mean = as.vector(X %*% beta), 
               sigma = solve(diag(N.obs) - C) %*% M))
  
  
  ################
  
  # phi.start <- runif(1, 1/ln, 1/l1)
  # var.start <- runif(1)
  # beta.start <- runif(N.beta)
  # m <- MakeADFun(
  #   data = list(W=W, X=X, Y=Y),
  #   parameters = list(beta=beta.start, phi=phi.start, variance=var.start),
  #   DLL = "carworks"
  # )
  # 
  
  m <- MakeADFun(
    data = list(W=W, X=X, Y=Y),
    parameters = list(beta=runif(2), sd=runif(1)),
    #DLL = "car_2"
    DLL = "car_simple"
  )
  #m$env$tracepar <- TRUE
  
  
  fit.new <- nlminb(m$par,m$fn,m$gr)
  #print(fit.new$par)
  #print(c(beta, phi, variance))
  
  # 
  # fit.new <- nlminb(m$par,m$fn,m$gr,
  #                   control=list(iter.max=300,eval.max=300))
  # report <- sdreport(m, getJointPrecision = 1)
  # 
  
  par.pred <- rbind(par.pred, fit.new$par)
  par.true <- rbind(par.true, c(beta, variance))
}

for(k in 1:3){
  plot(par.pred[, k], par.true[, k])
  abline(a=0, b=1)
  print(cor(par.pred[, k], par.true[, k]))
}
#print("phi")
#print(phi)

World2 <- fortify(admin2)
ids <- as.numeric(World2$id)
Y.pixel <- rep(Y - X%*%beta, sapply(0:21, function(i) sum(ids == i)))
World2$Y <- Y.pixel
p <- ggplot() +
  geom_polygon(data = World2, aes(x = long, y = lat, group = group, fill=Y.pixel),
               colour = "black", size = 0.5) + coord_fixed()
print(p)
