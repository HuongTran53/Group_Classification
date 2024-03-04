## --Filename: CF_antieigen.R
########################################
########################################
# Fixed covariance matrix 
p <- 5
fSigma <- diag(p); 
fSigma1 <- .75*(2/3 * diag(p) + 1/3 * rep(1, p) %*% t(rep(1, p)))
fSigma2 <- 2 * (.1 * diag(p) + .9 * rep(1, p) %*% t(rep(1, p)))
########################################
m <- 2e3
########################################
# 2 normal populations: Same Covariance matrix
n1 <- 200; n2 <- 200; p <- 5
mu1 <- rep(0, p); Sigma <- fSigma
mu2 <- rep(0, p) + 3; 
set.seed(1234)
X1 <- MASS::mvrnorm(n1, mu1, Sigma)
X2 <- MASS::mvrnorm(n2, mu2, Sigma)
# params <- list(mu1 = mu1, mu2 = mu2, Sigma = Sigma)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- MASS::mvrnorm(n1, mu1, Sigma)
  U2 <- MASS::mvrnorm(n2, mu2, Sigma)
  
  discr_type_1_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =1)$pop 
  discr_type_1_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =1)$pop 
  
  discr_type_2_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =2)$pop 
  discr_type_2_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =2)$pop 
  
  discr_type_3_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =3)$pop 
  discr_type_3_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =3)$pop 
  
  discr_type_4_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =4)$pop 
  discr_type_4_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =4)$pop 
  
  c(type1_11 = discr_type_1_1, 
    type1_22 = discr_type_1_2,
    type2_11 = discr_type_2_1, 
    type2_22 = discr_type_2_2,
    type3_11 = discr_type_3_1, 
    type3_22 = discr_type_3_2,
    type4_11 = discr_type_4_1, 
    type4_22 = discr_type_4_2
    )
})
re1 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
re1
########################################
########################################
# 2 normal populations: Same mean/Different Covariance matrix
n1 <- 200; n2 <- 200; p <- 5
mu1 <- rep(0, p); Sigma1 <- fSigma1
mu2 <- mu1; Sigma2 <- fSigma2
set.seed(1234)
X1 <- MASS::mvrnorm(n1, mu1, Sigma1)
X2 <- MASS::mvrnorm(n2, mu2, Sigma2)
# params <- list(mu1 = mu1, mu2 = mu2, Sigma = Sigma)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- MASS::mvrnorm(n1, mu1, Sigma1)
  U2 <- MASS::mvrnorm(n2, mu2, Sigma2)
  
  discr_type_1_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =1)$pop 
  discr_type_1_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =1)$pop 
  
  discr_type_2_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =2)$pop 
  discr_type_2_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =2)$pop 
  
  discr_type_3_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =3)$pop 
  discr_type_3_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =3)$pop 
  
  discr_type_4_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =4)$pop 
  discr_type_4_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =4)$pop 
  
  c(type1_11 = discr_type_1_1, 
    type1_22 = discr_type_1_2,
    type2_11 = discr_type_2_1, 
    type2_22 = discr_type_2_2,
    type3_11 = discr_type_3_1, 
    type3_22 = discr_type_3_2,
    type4_11 = discr_type_4_1, 
    type4_22 = discr_type_4_2
  )
})
re2 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
re2
########################################
########################################
# Kotz/t: Same Covariance matrix, different mean 
n1 <- 200; n2 <- 200; p <- 5
mu1 <- rep(0, p) + 3;  Sigma <- fSigma
mu2 <- rep(0, p); 
params <- list(mu1 = mu1, mu2 = mu2, Sigma = Sigma)
set.seed(1234)
source("~/OU /Dissertation/Normality Test/Kotz.R")
X1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma)
X2 <- mvtnorm::rmvt(n2, 5/7*Sigma, df = 7, delta = mu2)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma)
  U2 <- mvtnorm::rmvt(n2, 5/7*Sigma, df = 7, delta = mu2)
  
  discr_type_1_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =1)$pop 
  discr_type_1_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =1)$pop 
  
  discr_type_2_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =2)$pop 
  discr_type_2_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =2)$pop 
  
  discr_type_3_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =3)$pop 
  discr_type_3_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =3)$pop 
  
  discr_type_4_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =4)$pop 
  discr_type_4_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =4)$pop 
  
  c(type1_11 = discr_type_1_1, 
    type1_22 = discr_type_1_2,
    type2_11 = discr_type_2_1, 
    type2_22 = discr_type_2_2,
    type3_11 = discr_type_3_1, 
    type3_22 = discr_type_3_2,
    type4_11 = discr_type_4_1, 
    type4_22 = discr_type_4_2
  )
})
re3 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
re3
########################################
########################################
# Kotz/t: Same mean, difference covariance matrix
n1 <- 200; n2 <- 200; p <- 5
mu1 <- rep(0, p);  Sigma1 <- fSigma1
mu2 <- rep(0, p); Sigma2 <- fSigma2
set.seed(1234)
source("~/OU /Dissertation/Normality Test/Kotz.R")
# X1 <- MASS::mvrnorm(n1, mu1, Sigma1)
X1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma1)
X2 <- mvtnorm::rmvt(n2, 5/7*Sigma2, df = 7, delta = mu2)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma1)
  U2 <- mvtnorm::rmvt(n2, 5/7*Sigma2, df = 7, delta = mu2)
  
  discr_type_1_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =1)$pop 
  discr_type_1_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =1)$pop 
  
  discr_type_2_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =2)$pop 
  discr_type_2_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =2)$pop 
  
  discr_type_3_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =3)$pop 
  discr_type_3_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =3)$pop 
  
  discr_type_4_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =4)$pop 
  discr_type_4_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =4)$pop 
  
  c(type1_11 = discr_type_1_1,
    type1_22 = discr_type_1_2,
    type2_11 = discr_type_2_1, 
    type2_22 = discr_type_2_2,
    type3_11 = discr_type_3_1, 
    type3_22 = discr_type_3_2,
    type4_11 = discr_type_4_1, 
    type4_22 = discr_type_4_2
  )
})
re4 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean)
             )
re4
##########################################
##########################################
# A drawback of eigenstructure method 
# #Normal/Lomax: same mean/same covariance
# n1 <- 200; n2 <- 200; p <- 2
# a <- 3; theta = rep(1, p)
# S <- covLomax(a, theta)$Sigma
# mu <- covLomax(a, theta)$mu
# 
# set.seed(111)
# X1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
# X2 <- MASS::mvrnorm(n2, mu, S)
# 
# set.seed(999)
# set.seed(789)
# 
# re <- replicate(m, expr = {
#   
#   U1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
#   U2 <- MASS::mvrnorm(n2, mu, S)
#   
#   
#   discr_type_1_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =1)$pop 
#   discr_type_1_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =1)$pop 
#   
#   
#   discr_type_2_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =2)$pop 
#   discr_type_2_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =2)$pop 
#   
#   discr_type_3_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =3)$pop 
#   discr_type_3_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =3)$pop 
#   
#   discr_type_4_1 <- discr.antieigen(X1 = X1, X2 = X2, U = U1, type =4)$pop 
#   discr_type_4_2 <- discr.antieigen(X1 = X1, X2 = X2, U = U2, type =4)$pop 
#   
#   c(type1_11 = discr_type_1_1, 
#     type1_22 = discr_type_1_2,
#     type2_11 = discr_type_2_1, 
#     type2_22 = discr_type_2_2,
#     type3_11 = discr_type_3_1, 
#     type3_22 = discr_type_3_2,
#     type4_11 = discr_type_4_1, 
#     type4_22 = discr_type_4_2
#   )
# })
# 
# re4 <- rbind(apply(re == "Population 1", 1, mean), 
#              apply(re == "Can not decide", 1, mean), 
#              apply(re == "Population 2", 1, mean))
