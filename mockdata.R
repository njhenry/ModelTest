setwd("C://Users/scro3122/Documents/ModelTest")
#randomly set seed (so can reproduce any specific run)
seed <- runif(1)
set.seed(seed)

#make mock data where we:
#create rainfall and temp for each location
#make them affect cases on a two month timelag
#see if algorithm can detect this
cases_all <- c()
rain_all0 <- c()
rain_all1 <- c()
rain_all2 <- c()
temp_all0 <- c()
temp_all1 <- c()
temp_all2 <- c()
irs_all0 <- c()
irs_all1 <- c()
irs_all2 <- c()

N_locs <- 70

for(i in 1:N_locs){
  #just one location
  N_points <- 18
  l <- 0.1
  rain <- sin(seq(from=0, to=pi, length.out=N_points)) + l*rnorm(N_points)
  temp <- cos(seq(from=0, to=2*pi, length.out=N_points)) ^ 2 + l*rnorm(N_points)
  
  irs <- runif(N_points)
  l2 <- 0.5
  current_index <- 2 + (1:12)
  cases <- (rain + temp)[current_index - 2]^2 + rnorm(12, sd=0.5) - l2*irs[current_index - 1]
  rain0 <- rain[current_index]
  rain1 <- rain[current_index - 1]
  rain2 <- rain[current_index - 2]
  temp0 <- temp[current_index]
  temp1 <- temp[current_index - 1]
  temp2 <- temp[current_index - 2]
  irs0 <- irs[current_index]
  irs1 <- irs[current_index - 1]
  irs2 <- irs[current_index - 2]
  
  cases_all <- c(cases_all, cases)
  rain_all0 <- c(rain_all0, rain0)
  rain_all1 <- c(rain_all1, rain1)
  rain_all2 <- c(rain_all2, rain2)
  temp_all0 <- c(temp_all0, temp0)
  temp_all1 <- c(temp_all1, temp1)
  temp_all2 <- c(temp_all2, temp2)
  irs_all0 <- c(irs_all0, irs0)
  irs_all1 <- c(irs_all1, irs1)
  irs_all2 <- c(irs_all2, irs2)
}
obsDat <- list(cases_all, 
               rain_all0, rain_all1, rain_all2, 
               temp_all0, temp_all1, temp_all2, 
               irs_all0, irs_all1, irs_all2)
names(obsDat) <- c("cases",
                   paste0("rain-", 0:2),
                   paste0("temp-", 0:2),
                   paste0("IRS-", 0:2))
n_vars <- length(obsDat)
G_0 <- matrix(1, nrow=n_vars, ncol=n_vars) - diag(n_vars)
G_0[1, ] <- 0

source("runpc.R")
pc <- runpc(obsDat, alpha=0.7, pw=FALSE, G_0=G_0)
plot.minimal.n(pc[[1]], 1, names(obsDat), 5)
