## --Filename: CF_MaxSep.R
## Data experiment for majority class rule
########################################
# Fixed covariance matrix 
p <- 5
fSigma <- diag(p); 
fSigma1 <- .75*(2/3 * diag(p) + 1/3 * rep(1, p) %*% t(rep(1, p)))
fSigma2 <- 2 * (.1 * diag(p) + .9 * rep(1, p) %*% t(rep(1, p)))
########################################
m <- 2e3
# 2 normal populations: Same Covariance matrix
n1 <- 200; n2 <- 200; p <- 5
mu1 <- rep(0, p); Sigma <- fSigma
mu2 <- rep(0, p) + 3; 
set.seed(1234)
X1 <- MASS::mvrnorm(n1, mu1, Sigma)
X2 <- MASS::mvrnorm(n2, mu2, Sigma)
params <- list(mu1 = mu1, mu2 = mu2, Sigma = Sigma)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- MASS::mvrnorm(n1, mu1, Sigma)
  U2 <- MASS::mvrnorm(n2, mu2, Sigma)

  re.linear_known1 <- discr.linear(U = U1, params = params)$pop 
  re.linear_known2 <- discr.linear(U = U2, params = params)$pop 

  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 

  re.MaxSep_known1 <- discr.MaxSep(U = U1,params = params)$pop 
  re.MaxSep_known2 <- discr.MaxSep(U = U2, params = params)$pop 

  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop

  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop

  tmp <- c(linear_known11 = re.linear_known1,
            linear_known22 = re.linear_known2,
            linear11 = re.linear1,
            linear22 = re.linear2,
            MaxSep_known11 = re.MaxSep_known1,
            MaxSep_known22 = re.MaxSep_known2,
            MaxSep11 = re.MaxSep1,
            MaxSep22 = re.MaxSep2,
            avgD11 = re.avgD1, 
            avgD22 = re.avgD2, 
            leverage11 = re.leverage1, 
            leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
})

re1 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re1)
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
  
  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 
  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
  
  tmp <- c(
    linear11 = re.linear1,
    linear22 = re.linear2,
    MaxSep11 = re.MaxSep1,
    MaxSep22 = re.MaxSep2,
    avgD11 = re.avgD1, 
    avgD22 = re.avgD2, 
    leverage11 = re.leverage1, 
    leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
})

re2 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re2)
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
  
  re.linear_known1 <- discr.linear(U = U1, params = params)$pop 
  re.linear_known2 <- discr.linear(U = U2, params = params)$pop 
  
  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 
  
  re.MaxSep_known1 <- discr.MaxSep(U = U1,params = params)$pop 
  re.MaxSep_known2 <- discr.MaxSep(U = U2, params = params)$pop 
  
  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
  
  tmp <- c(linear_known11 = re.linear_known1,
           linear_known22 = re.linear_known2,
           linear11 = re.linear1,
           linear22 = re.linear2,
           MaxSep_known11 = re.MaxSep_known1,
           MaxSep_known22 = re.MaxSep_known2,
           MaxSep11 = re.MaxSep1,
           MaxSep22 = re.MaxSep2,
           avgD11 = re.avgD1, 
           avgD22 = re.avgD2, 
           leverage11 = re.leverage1, 
           leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
  
})

re3 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re3)
########################################
########################################
# Kotz/t: Same mean, different covariance matrix
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
  
  
  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 
  
  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
  
  tmp <- c(
    # linear_known11 = re.linear_known1,
    # linear_known22 = re.linear_known2,
    linear11 = re.linear1,
    linear22 = re.linear2,
    # MaxSep_known11 = re.MaxSep_known1,
    # MaxSep_known22 = re.MaxSep_known2,
    MaxSep11 = re.MaxSep1,
    MaxSep22 = re.MaxSep2,
    avgD11 = re.avgD1, 
    avgD22 = re.avgD2, 
    leverage11 = re.leverage1, 
    leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
  
})

re4 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re4)

#######################################
########################################
# Lomax/t: Same mean, diff covariance matrix
n1 <- 200; n2 <- 200; p <- 5
a <- 3; theta = rep(1, p)
S1 <- covLomax(a, theta)$Sigma
mu <- covLomax(a, theta)$mu
S2 <- fSigma2
set.seed(111)
source("~/OU /Dissertation/Normality Test/Kotz.R")
# X1 <- MASS::mvrnorm(n1, mu1, Sigma1)
X1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
X2 <- mvtnorm::rmvt(n2, 5/7*S2, df = 7, delta = mu)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
  U2 <- mvtnorm::rmvt(n2, 5/7*S2, df = 7, delta = mu)
  
  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 
  
  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
  
  tmp <- c(
    linear11 = re.linear1,
    linear22 = re.linear2,
    MaxSep11 = re.MaxSep1,
    MaxSep22 = re.MaxSep2,
    avgD11 = re.avgD1, 
    avgD22 = re.avgD2, 
    leverage11 = re.leverage1, 
    leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
  
})

re5 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re5)
###################################################
# Lomax/t: Same Covariance matrix, different mean 
n1 <- 200; n2 <- 200; p <- 5
a <- 3; theta = rep(1, p)
S <- covLomax(a, theta)$Sigma
mu1 <- covLomax(a, theta)$mu
mu2 <- rep(3, p)

set.seed(111)
source("~/OU /Dissertation/Normality Test/Kotz.R")
# X1 <- MASS::mvrnorm(n1, mu1, Sigma1)
X1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
X2 <- mvtnorm::rmvt(n2, 5/7*S, df = 7, delta = mu2)
params <- list(mu1 = mu1, mu2 = mu2, Sigma = S)
set.seed(1234)
source("~/OU /Dissertation/Normality Test/Kotz.R")
set.seed(789)
re <- replicate(m, expr = {
  
  U1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
  U2 <- mvtnorm::rmvt(n2, 5/7*S, df = 7, delta = mu2)
  
  re.linear_known1 <- discr.linear(U = U1, params = params)$pop 
  re.linear_known2 <- discr.linear(U = U2, params = params)$pop 
  
  re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2, true.params = F)$pop 
  
  re.MaxSep_known1 <- discr.MaxSep(U = U1,params = params)$pop 
  re.MaxSep_known2 <- discr.MaxSep(U = U2, params = params)$pop 
  
  re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
  
  re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
  
  tmp <- c(linear_known11 = re.linear_known1,
           linear_known22 = re.linear_known2,
           linear11 = re.linear1,
           linear22 = re.linear2,
           MaxSep_known11 = re.MaxSep_known1,
           MaxSep_known22 = re.MaxSep_known2,
           MaxSep11 = re.MaxSep1,
           MaxSep22 = re.MaxSep2,
           avgD11 = re.avgD1, 
           avgD22 = re.avgD2, 
           leverage11 = re.leverage1, 
           leverage22 = re.leverage2)
  # Randomly assign to each population if Can not decide.
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
  
})
re6 <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Population 2", 1, mean))
print(re6)
# ##############################################
# ##############################################
# # Do Not run 
# # Normal/Lomax: same covariance
# n1 <- 200; n2 <- 200; p <- 5
# a <- 3; theta = rep(1, p)
# S <- covLomax(a, theta)$Sigma
# mu <- covLomax(a, theta)$mu
# 
# set.seed(111)
# X1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
# X2 <- MASS::mvrnorm(n2, mu, S)
# set.seed(999)
# re <- replicate(m, expr = {
#   
#   U1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
#   U2 <- MASS::mvrnorm(n2, mu, S)
# 
# 
#   re.linear_known1 <- discr.linear(U = U1,
#                                    mu1 = mu, mu2 = mu, Sigma = S)$pop
#   re.linear_known2 <- discr.linear(U = U2,
#                                    mu1 = mu, mu2 =mu, Sigma = S)$pop
# 
#   re.linear1 <- discr.linear(X1 = X1, X2 = X2, U = U1)$pop
#   re.linear2 <- discr.linear(X1 = X1, X2 = X2, U = U2)$pop
# 
# 
#   re.MaxSep_known1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1,
#                                    mu1 = mu, mu2 =mu, Sigma = S)$pop
#   re.MaxSep_known2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2,
#                                    mu1 = mu, mu2 =mu, Sigma = S)$pop
# 
#   re.MaxSep1 <- discr.MaxSep(X1 = X1, X2 = X2, U = U1)$pop
#   re.MaxSep2 <- discr.MaxSep(X1 = X1, X2 = X2, U = U2)$pop
# 
#   re.avgD1 <- discr.averageD(X1 = X1, X2 = X2, U = U1)$pop 
#   re.avgD2 <- discr.averageD(X1 = X1, X2 = X2, U = U2)$pop 
#   
#   re.leverage1 <- discr.leverage(X1 = X1, X2 = X2, U = U1)$pop
#   re.leverage2 <- discr.leverage(X1 = X1, X2 = X2, U = U2)$pop
#   
#   c(
#     linear_known11 = re.linear_known1,
#     linear_known22 = re.linear_known2,
#     linear11 = re.linear1,
#     linear22 = re.linear2,
#     MaxSep_known11 = re.MaxSep_known1,
#     MaxSep_known22 = re.MaxSep_known2,
#     MaxSep11 = re.MaxSep1,
#     MaxSep22 = re.MaxSep2,
#     avgD11 = re.avgD1, 
#     avgD22 = re.avgD2, 
#     leverage11 = re.leverage1, 
#     leverage22 = re.leverage2)
# })
# 
# re4 <- rbind(apply(re == "Population 1", 1, mean), 
#              apply(re == "Can not decide", 1, mean), 
#              apply(re == "Population 2", 1, mean))
# 
# 
# miss_points_normal_normal <- re1
# save(miss_points_normal_normal, file = "miss_points_normal_normal.RDATA")
# miss_points_kotz_t <- re2
# save(miss_points_kotz_t, file = "miss_points_kotz_t.RDATA")
# miss_points_kotz_t2 <- re3
# save(miss_points_kotz_t2, file = "miss_points_kotz_t2.RDATA")
# miss_points_lomax_normal <- re4
# save(miss_points_lomax_normal, file = "miss_points_lomax_normal.RDATA")
















