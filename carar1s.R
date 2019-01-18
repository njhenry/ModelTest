setwd("C:/Users/scro3122/Documents/ModelTest")
library(rgeos)
library(raster)
library(mvtnorm)
library(rgdal)
library(ggplot2)
library(TMB)


admin2 <- shapefile("C:/Users/scro3122/Documents/Madagascar_Admin_shapefiles/Malareo_Region.shp")
W <- gTouches(admin2, byid=TRUE) * 1

N.reps <- 100
beta.true <- c()
beta.pred <- c()
beta.upper <- c()
beta.lower <- c()
beta.upper.mh <- c()
beta.lower.mh <- c()


compile("carar1.cpp")
dyn.load(dynlib("carar1"))


N.beta <- 4
for(k in 1:N.reps){

  print(k)
  
  #######create fake data
  N.obs <- dim(W)[1]
  phi <- runif(1)
  #phi <- 1

  N.times <- 24
  X <- matrix(runif(N.beta * N.obs), ncol=N.beta)
  beta <- runif(N.beta, -1, 1)
  C <- W / N.obs
  variance <- runif(1)
  M <- diag(N.obs) * variance
  #create first time obs
  #Y.int <- as.vector(rmvnorm(1, mean = as.vector(X %*% beta), sigma = solve(diag(N.obs) - C) %*% M))
  Y.int <- as.vector(rmvnorm(1, mean = as.vector(X %*% beta), sigma = solve(diag(N.obs) - C) %*% M))
  Y <- Y.int
  Y.old <- Y
  
  for(i in 2:N.times){
    #Y.new <- as.vector(rmvnorm(1, mean = phi * Y.old + as.vector(X %*% beta), sigma = solve(diag(N.obs) - C) %*% M))
    Y.new <- as.vector(rmvnorm(1, mean = phi * (Y.old - as.vector(X %*% beta))  + as.vector(X %*% beta), 
                               sigma = solve(diag(N.obs) - C) %*% M))
    #Y.new <- as.vector(phi * (Y.old - as.vector(X %*% beta))  + as.vector(X %*% beta))
    Y <- rbind(Y, Y.new)
    Y.old <- Y.new
  }
  
  

  
  
  
  # ######fit
  m <- MakeADFun(
    data = list(sd_phi=0.1, W=W, X=X, Y=Y),
    parameters = list(beta=runif(N.beta), sd=runif(1), inv_phi=runif(1)),
    #DLL = "car_2",
    DLL = "carar1"
  )
  #m$env$tracepar <- TRUE
  

  
  fit.new <- nlminb(m$par,m$fn,m$gr)
  
  # ##MH
  # N.par <- length(fit.new$par)
  # N.mh <- 1000
  # par.samples <- matrix(NA, nrow=N.mh, ncol=N.par)
  # par.samples[1, ] <- fit.new$par
  # k <- 0.005
  # for(i in 2:N.mh){
  #   par.current <- par.samples[i-1, ]
  #   par.candidate <- par.current + k*rnorm(N.par)
  # 
  #   l <- m$fn(par.candidate) / m$fn(par.current)
  #   #print(m$fn(par.candidate))
  #   #print(l)
  # 
  #   if(l > runif(1)){
  #     par.samples[i, ] <- par.candidate
  #   }else{
  #     par.samples[i, ] <- par.current
  #   }
  # }
  # mh.sd <- apply(par.samples, 2, sd)[1:N.beta]
  
  
  #try and compute standard deviations manually
  sd.beta <- sqrt(diag(solve(m$he(fit.new$par)))[1:N.beta])
  
  #sds <- sdreport(m, getJointPrecision = 0)
  # sd.beta <- sqrt(diag(sds$cov.fixed)[1:N.beta])
  
  alpha <- 0.8
  
  # beta.pred.upper.mh <- fit.new$par[1:N.beta] + qnorm((1+alpha)/2) * mh.sd
  # beta.pred.lower.mh <- fit.new$par[1:N.beta] - qnorm((1+alpha)/2) * mh.sd
  # beta.upper.mh <- rbind(beta.upper.mh, beta.pred.upper.mh)
  # beta.lower.mh <- rbind(beta.lower.mh, beta.pred.lower.mh)
 
  beta.pred.upper <- fit.new$par[1:N.beta] + qnorm((1+alpha)/2) * sd.beta
  beta.pred.lower <- fit.new$par[1:N.beta] - qnorm((1+alpha)/2) * sd.beta
  
  
  beta.true <- rbind(beta.true, beta)
  beta.pred <- rbind(beta.pred, fit.new$par[1:N.beta])
  beta.upper <- rbind(beta.upper, beta.pred.upper)
  beta.lower <- rbind(beta.lower, beta.pred.lower)
}

for(k in 1:N.beta){
  b.df <- data.frame(true=beta.true[, k], pred=beta.pred[, k],
                     upper=beta.upper[, k], lower=beta.lower[, k], iter=1:length(beta.true[, k]),
                     contained = (beta.lower[, k] < beta.true[, k]) & (beta.upper[, k] > beta.true[, k]),
                     contzero = (beta.lower[, k] * beta.upper[, k] < 0))
  n <- length(beta.pred[, k])
  # b2.df <- data.frame(vals=c(beta.true[, k], beta.pred[, k], beta.lower[, k], beta.upper[, k]),
  #                     iter=1:n,
  #                     type=rep(c("true", "pred", "lower", "upper"), each = n))
  
  b2.df <- data.frame(vals=c(beta.true[, k], beta.lower[, k], beta.upper[, k]),
                      iter=1:n,
                      type=rep(c("true", "lower", "upper"), each = n))
  
  b3.df <- b.df[order(b.df$true), ]
  b3.df$iter <- 1:n
  
  # p <- ggplot(b2.df) + geom_point(aes(x=iter, y=vals, color=type))
  # print(p)
  
  p <- ggplot(b3.df) + geom_point(aes(x=iter, y=true, color=contained)) + 
    geom_errorbar(aes(x=iter, ymin=lower, ymax=upper), alpha = 0.3) 
  print(p)
  
  print(sum((beta.lower[, k] < beta.true[, k]) & (beta.upper[, k] > beta.true[, k])) / n)
}

# 
# #mh
# for(k in 1:2){
#   b.df <- data.frame(true=beta.true[, k], pred=beta.pred[, k],
#                      upper=beta.upper.mh[, k], lower=beta.lower.mh[, k], iter=1:length(beta.true[, k]),
#                      contained = (beta.lower.mh[, k] < beta.true[, k]) & (beta.upper.mh[, k] > beta.true[, k]),
#                      contzero = (beta.lower.mh[, k] * beta.upper.mh[, k] < 0))
#   n <- length(beta.pred[, k])
#   # b2.df <- data.frame(vals=c(beta.true[, k], beta.pred[, k], beta.lower[, k], beta.upper[, k]),
#   #                     iter=1:n,
#   #                     type=rep(c("true", "pred", "lower", "upper"), each = n))
#   
#   # b2.df <- data.frame(vals=c(beta.true[, k], beta.lower[, k], beta.upper[, k]),
#   #                     iter=1:n,
#   #                     type=rep(c("true", "lower", "upper"), each = n))
#   
#   b3.df <- b.df[order(b.df$true), ]
#   b3.df$iter <- 1:n
#   
#   # p <- ggplot(b2.df) + geom_point(aes(x=iter, y=vals, color=type))
#   # print(p)
#   
#   p <- ggplot(b3.df) + geom_point(aes(x=iter, y=true, color=contained)) + 
#     geom_errorbar(aes(x=iter, ymin=lower, ymax=upper), alpha = 0.3) 
#   print(p)
#   
#   print(sum((beta.lower.mh[, k] < beta.true[, k]) & (beta.upper.mh[, k] > beta.true[, k])) / n)
# }



# for(k in 1:(dim(beta.true)[2])){
#   plot(beta.true[, k], beta.pred[, k])
#   abline(a=0, b=1)
#   print(cor(beta.true[, k], beta.pred[, k]))
# }

# #visualize a few times
# World2 <- fortify(admin2)
# ids <- as.numeric(World2$id)
# 
# for(i in 1:3){
#   #Y.pixel <- rep(Y[i, ] - X%*%beta, sapply(0:21, function(i) sum(ids == i)))
#   Y.pixel <- rep(Y[i, ], sapply(0:21, function(i) sum(ids == i)))
#   World2$Y <- Y.pixel
#   p <- ggplot() +
#     geom_polygon(data = World2, aes(x = long, y = lat, group = group, fill=Y.pixel),
#                  colour = "black", size = 0.5) + coord_fixed() + ggtitle(paste0(i))
#   print(p)
# }
# 
