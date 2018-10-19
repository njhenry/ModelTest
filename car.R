setwd("C:/Users/scro3122/Documents/ModelTest/")
N.obs <- 100
N.beta <- 2
beta <- runif(N.beta)
X <- matrix(runif(N.obs * N.beta), nrow = N.obs)

Y <- as.vector(X %*% beta) + rnorm(N.obs, sd=0.1)

W <- matrix(NA, nrow=N.obs, ncol=N.obs)

for(i in 1:N.obs){
  for(j in 1:i){
    W[i, j] <- sample(c(0, 1), 1)
    W[j, i] <- W[i, j]
  }
}

# W <- matrix(c(0,1,1,0,
#               1,0,0,1,
#               1,0,0,1,
#               0,1,1,0), nrow = 4, byrow=TRUE)


library(TMB)


compile("car.cpp")
dyn.load(dynlib("car"))

m <- MakeADFun(
  data = list(Y=runif(N.obs), X=X, Y=Y, W=W),
  parameters = list(sigma=runif(1), phi=runif(1), beta=runif(N.beta)),
  DLL = "car"
)

fit.new <- nlminb(m$par,m$fn,m$gr)