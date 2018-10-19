library(Rcpp)
library(FSIC)
#sourceCpp("C:/Users/scro3122/Documents/ModelTest/fsicbootstrap.cpp")


# fstest_b <- function(x, y, alpha = 0.05, test = FALSE, J = 20, nBt = 200){
#   if(length(x) != length(y)) stop("vectors are of different lengths")
# 
#   #remove entries where either variable is NA
#   keep <- !(is.na(x)|is.na(y))
#   x <- x[keep]
#   y <- y[keep]
# 
#   N <- length(x)
# 
#   v <- rnorm(J)
#   w <- rnorm(J)
#   sx <- 1
#   sy <- 1
#   pval <- fstest_debug(x, y, v, w, nBt)
#   if(test) return(pval < alpha)
#   return(pval)
# }

N <- 1000
x <- runif(N)
y <- runif(N)
z <- x + rnorm(N)

# print(fstest_b(x, y))
# print(fstest_b(x, x))
# print(fstest_b(x, z))

x2 <- rep(0, N)
x2[2] <- 1
y2 <- rep(0, N)
y2[5] <- 1

print(fstest(x2, y2))

# N <- 1000
# x <- runif(N)
# y <- runif(N)
# z <- x + rnorm(N)
# 
# print(fstest(x, y))
# print(fstest(x, x))
# print(fstest(x, z))