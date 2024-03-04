## -- Filename: CF_NP.R
m <- 2e3; 
n1 <- 200; n2 <- 200; p <- 5
source("~/OU /Dissertation/Normality Test/Kotz.R")
######################################
######################################
# 2 normal populations: Same Covariance matrix
mu1 <- rep(0, p); Sigma <- fSigma
mu2 <- rep(0, p) + 3; 
f1 <- function(u) mvtnorm::dmvnorm(u, mean = mu1, sigma = Sigma)
f2 <- function(u) mvtnorm::dmvnorm(u, mean = mu2, sigma = Sigma)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- MASS::mvrnorm(n1, mu1, Sigma)
  U2 <- MASS::mvrnorm(n2, mu2, Sigma)
  re1.known <- ifelse(
    NP(U1, lst_f = c(f1, f2)) == 1, 
    "Population 1", "Population 2"
  )
  re2.known <- ifelse(
    NP(U2, lst_f = c(f1, f2)) == 1, 
    "Population 1", "Population 2"
  )
  c(NP_known11 = re1.known,
    NP_known22 = re2.known
  )
})
re1_NP <- rbind(apply(re == "Population 1", 1, mean), 
                apply(re == "Population 2", 1, mean))
print(re1_NP)
######################################
######################################
# 2 normal populations: Same mean/Different Covariance matrix
mu1 <- rep(0, p); Sigma1 <- fSigma1
mu2 <- mu1; Sigma2 <- fSigma2
set.seed(789)
f1 <- function(u) mvtnorm::dmvnorm(u, mean = mu1, sigma = Sigma1)
f2 <- function(u) mvtnorm::dmvnorm(u, mean = mu2, sigma = Sigma2)

set.seed(789)
re <- replicate(m, expr = {
  U1 <- MASS::mvrnorm(n1, mu1, Sigma1)
  U2 <- MASS::mvrnorm(n2, mu2, Sigma2)
  re1.known <- ifelse(
    NP(U1, lst_f = c(f1, f2)) == 1, 
    "Population 1", "Population 2"
  )
  re2.known <- ifelse(
    NP(U2, lst_f = c(f1, f2)) == 1, 
    "Population 1", "Population 2"
  )
  c(NP_known11 = re1.known,
    NP_known22 = re2.known
  )
})
re2_NP <- rbind(apply(re == "Population 1", 1, mean), 
                apply(re == "Population 2", 1, mean))
print(re2_NP)
########################################
########################################
# Kotz/t: Same Covariance matrix, different mean
mu1 <- rep(0, p) + 3;  Sigma <- fSigma
mu2 <- rep(0, p); 
f1 <- function(u) log(dmvKotz(u, mu = mu1, Sigma = 1/(p+1)*Sigma))
f2 <- function(u) mvtnorm::dmvt(u, 5/7*Sigma, df = 7, delta = mu2, log = T)
set.seed(789)
re <- replicate(m, expr = {
  U1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma)
  U2 <- mvtnorm::rmvt(n2, 5/7*Sigma, df = 7, delta = mu2)
  re1.known <- ifelse(
    NP(U1, lst_f = c(f1, f2), use.log = F) == 1, 
    "Population 1", "Population 2"
  )
  re2.known <- ifelse(
    NP(U2, lst_f = c(f1, f2), use.log = F) == 1,
    "Population 1", "Population 2"
  )
  # For U1:
  re1 <- mleKotz(colMeans(U1), cov(U1), U1)
  hf1 <- function(u) log(dmvKotz(u, mu = re1$mu, Sigma = re1$Sigma))
  re2 <- Compositional::multivt(U1)
  hf2 <- function(u) mvtnorm::dmvt(u, 
                                   ((re2$df-2)/re2$df * re2$covariance),
                                   df = re2$df, delta = re2$center,log = T
  )
  re1_unknown <- ifelse(
    NP(U1, lst_f = c(hf1, hf2), use.log = F) == 1,
    "Population 1", "Population 2"
  )
  # For U2:
  re1 <- mleKotz(colMeans(U2), cov(U2), U2)
  hf1 <- function(u) log(dmvKotz(u, mu = re1$mu, Sigma = re1$Sigma))
  re2 <- Compositional::multivt(U2)
  hf2 <- function(u) mvtnorm::dmvt(u,
                                   ((re2$df-2)/re2$df * re2$covariance),
                                   df = re2$df, delta = re2$center, log = T
  )
  re2_unknown <- ifelse(NP(U2, lst_f = c(hf1, hf2), use.log = F) == 1, 
                        "Population 1", "Population 2"
  )
  c(NP_known11 = re1.known,
    NP_known22 = re2.known,
    NP11 = re1_unknown,
    NP22 = re2_unknown
  )
})
re3_nP <- rbind(apply(re == "Population 1", 1, mean), 
                apply(re == "Population 2", 1, mean))
print(re3_nP)
########################################
########################################
# Kotz/t: Same mean, difference covariance matrix
mu1 <- rep(0, p);  Sigma1 <- fSigma1
mu2 <- rep(0, p); Sigma2 <- fSigma2
# source("~/OU /Dissertation/Normality Test/Kotz.R")
f1 <- function(u) log(dmvKotz(u, mu = mu1, Sigma = 1/(p+1)*Sigma1))
f2 <- function(u) mvtnorm::dmvt(u, 5/7*Sigma2, df = 7, delta = mu2, log = T)
set.seed(789)
re <- replicate(m, expr = {
  
  U1 <- rmvKotz(n = n1, mu = mu1, Sigma = 1/(p+1)*Sigma1)
  U2 <- mvtnorm::rmvt(n2, 5/7*Sigma2, df = 7, delta = mu2)
  
  re1.known <- ifelse(
    NP(U1, lst_f = c(f1, f2), use.log = F) == 1,
    "Population 1", "Population 2"
  )
  re2.known <- ifelse(
    NP(U2, lst_f = c(f1, f2), use.log = F) == 1, 
    "Population 1", "Population 2"
  )
  # For U1:
  re1 <- mleKotz(colMeans(U1), cov(U1), U1)
  hf1 <- function(u) log(dmvKotz(u, mu = re1$mu, Sigma = re1$Sigma))
  re2 <- Compositional::multivt(U1)
  hf2 <- function(u) mvtnorm::dmvt(u, 
                                   ((re2$df-2)/re2$df * re2$covariance), 
                                   df = re2$df, delta = re2$center, log = T
  )
  re1_unknown <- ifelse(
    NP(U1, lst_f = c(hf1, hf2), use.log = F) == 1, 
    "Population 1", "Population 2"
  )
  # For U2:
  re1 <- mleKotz(colMeans(U2), cov(U2), U2)
  hf1 <- function(u) log(dmvKotz(u, mu = re1$mu, Sigma = re1$Sigma))
  re2 <- Compositional::multivt(U2)
  hf2 <- function(u) mvtnorm::dmvt(u,
                                   ((re2$df-2)/re2$df * re2$covariance), 
                                   df = re2$df, delta = re2$center, log = T
  )
  re2_unknown <- ifelse(
    NP(U2, lst_f = c(hf1, hf2), use.log = F) == 1, 
    "Population 1", "Population 2"
  )
  c(NP_known11 = re1.known,
    NP_known22 = re2.known,
    NP11 = re1_unknown,
    NP22 = re2_unknown
  )
})

re4_NP <- rbind(apply(re == "Population 1", 1, mean), 
                apply(re == "Population 2", 1, mean))
print(re4_NP)
##############################################
##############################################
# A drawback of eigenstructure method 
# # For lomax:
# mleLomax <- function(parm1, parm2, X){
#   n <- nrow(X); p <- ncol(X)
#   params0 <- c(parm1, parm2)
#   loglik.lomax <- function(params) {
#     ll <- sum(NonNorMvtDist::dmvlomax(X, 
#                                       parm1 = params[1], parm2 = params[-1],
#                                       log = TRUE
#     )
#     )
#     return(ll)
#   }
#   
#   est = constrOptim(
#     params0, f = loglik.lomax, grad = NULL, 
#     ui = diag(p + 1), ci = rep(0, p + 1), control = list(fnscale = -1)
#   )
#   
#   is.converge <- ifelse(est$convergence == 0, T, F)
#   return(list(
#     parm1 = est$par[1], parm2 = est$par[-1], Is.converge = is.converge)
#   )
# }
# 
# # Normal/Lomax: same covariance
# n1 <- 200; n2 <- 200; p <- 5
# a <- 3; theta = rep(1, p)
# S <- covLomax(a, theta)$Sigma
# mu <- covLomax(a, theta)$mu
# 
# set.seed(111)
# f1 <- function(u) NonNorMvtDist::dmvlomax(u, parm1 = a, parm2 = theta)
# f2 <- function(u) mvtnorm::dmvnorm(u, mean = mu, sigma = S)
# # m <- 100
# set.seed(999)
# re <- replicate(m, expr = {
#   
#   U1 <- NonNorMvtDist::rmvlomax(n1, parm1 = a, parm2 = theta)
#   U2 <- MASS::mvrnorm(n2, mu, S)
#   
#   re1.known <- ifelse(
#     NP(U1, lst_f = c(f1, f2)) == 1, "Population 1", "Population 2"
#   )
#   re2.known <- ifelse(
#     NP(U2, lst_f = c(f1, f2)) == 1, "Population 1", "Population 2"
#   )
#   # For U1:
#   re1 <- mleLomax(parm1 = 1, parm2 = rep(1, p), U1)
#   hf1 <- function(u) NonNorMvtDist::dmvlomax(u, 
#                                              parm1 = re1$parm1,
#                                              parm2 = re1$parm2
#   )
#   hf2 <- function(u) mvtnorm::dmvnorm(u, mean = colMeans(U1), sigma = cov(U1))
#   re1_unknown <- ifelse(
#     NP(U1, lst_f = c(hf1, hf2), use.log = F) == 1, 
#     "Population 1", "Population 2"
#   )
#   
#   # For U2:
#   re2 <- mleLomax(parm1 = 1, parm2 = rep(1, p), U2)
#   hf1 <- function(u) NonNorMvtDist::dmvlomax(u, 
#                                              parm1 = re2$parm1,
#                                              parm2 = re2$parm2
#   )
#   hf2 <- function(u) mvtnorm::dmvnorm(u, mean = colMeans(U2), sigma = cov(U2))
#   re2_unknown <- ifelse(
#     NP(U2, lst_f = c(hf1, hf2), use.log = F) == 1,
#     "Population 1", "Population 2"
#   )
#   c(NP_known11 = re1.known,
#     NP_known22 = re2.known,
#     NP11 = re1_unknown,
#     NP22 = re2_unknown
#   )
#   
# })
# 
# re4 <- rbind(apply(re == "Population 1", 1, mean),
#              apply(re == "Population 2", 1, mean))
# print(re4)
