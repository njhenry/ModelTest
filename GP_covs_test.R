setwd("C:/Users/scro3122/Documents/ModelTest")
library(INLA)
library(TMB)
library(RandomFields)
library(ggplot2)

M <- 1000
coord <- cbind(runif(M), runif(M))

N <- 950
hold.out <- sample.int(M, N)

a <- RMmatern(0.5)
cov <- RFsimulate(a, x=coord[,1], y=coord[,2])

a2 <- RMmatern(0.5)
coeff <- RFsimulate(a2, x=coord[,1], y=coord[,2])

response <- coeff$variable1 * cov$variable1 + rnorm(M, sd = 0.01)

df <- data.frame(x=coord[,1], y=coord[,2], cov=cov$variable1, coeff=coeff$variable1, response=response)

p <- ggplot(df, aes(x=x, y=y, color=cov)) + geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10))
print(p)

p <- ggplot(df, aes(x=x, y=y, color=coeff)) + geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10))
print(p)

p <- ggplot(df, aes(x=x, y=y, color=response)) + geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10))
print(p)

mesh <- inla.mesh.2d(loc = coord, cutoff = 0.1, max.n = 500) 
plot(mesh)
spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]
A <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coord[-hold.out, ]))
A_hold <- inla.spde.make.A(mesh=mesh, loc=as.matrix(coord[hold.out, ]))
n_s <- nrow(spde$M0)


compile("GP_covs_test.cpp")
dyn.load(dynlib("GP_covs_test"))



m <- MakeADFun(
  data = list(response=response[-hold.out], spde=spde, cov=cov$variable1[-hold.out], A=A),
  parameters = list(gp=runif(n_s, 0, 10), log_kappa=2.5, log_tau=0.0,
                    sigma=1.00),
  random = "gp",
  DLL = "GP_covs_test"
)



fit.new <- nlminb(m$par,m$fn,m$gr,
                  control=list(iter.max=300,eval.max=300))
print(fit.new)



rep <- sdreport(m, getJointPrecision = 1)


field.out <- rep$par.random
gp <- A_hold %*% field.out

pred <- as.vector(gp) * cov$variable1[hold.out]

plot(pred, response[hold.out])


df.used <- data.frame(x=coord[-hold.out,1], y=coord[-hold.out,2], 
                      cov=cov$variable1[-hold.out], coeff=coeff$variable1[-hold.out], response=response[-hold.out])

p <- ggplot(df.used, aes(x=x, y=y, color=response)) + geom_point() +
  scale_colour_gradientn(colours = terrain.colors(10))
print(p)


# 
# df <- data.frame(x=coord[,1], y=coord[,2], cov=cov$variable1, coeff=coeff$variable1, 
#                  response=response, pred_field=as.vector(gp))
# 
# p <- ggplot(df, aes(x=x, y=y, color=pred_field)) + geom_point() +
#   scale_colour_gradientn(colours = terrain.colors(10))
# print(p)
# 
# 
# p <- ggplot(df, aes(x=x, y=y, color=coeff)) + geom_point() +
#   scale_colour_gradientn(colours = terrain.colors(10))
# print(p)

