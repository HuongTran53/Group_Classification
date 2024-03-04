## --Filename: covLomax.R
covLomax <- function(a, theta){
  # Get covariance matrix of multivariate lomax distribution
  # a >= 2 to ensure finite covariance 
  p <- length(theta)
  S <- array(dim = c(p, p))
  mu <- c()
  for (i in 1:p){
    tmpi <- (gamma(a - 1)*gamma(2))/(gamma(a) * theta[i])
    mu <- c(mu, tmpi)
    tmpi2 <- (gamma(a - 2)*gamma(3))/(gamma(a)*theta[i]^2)
    S[i, i] <- tmpi2 - tmpi*tmpi
    
    for (j in (i+1):p){
      if (j > p){break}
      tmpj <- (gamma(a - 1)*gamma(2))/(gamma(a) * theta[j])
      tmp1 <- (gamma(a - 2)*gamma(2)^2)/(gamma(a) * theta[i] *theta[j])
      S[i, j] <- S[j, i] <- tmp1 - tmpi*tmpj
    }
  }
  
  return(list(mu = mu, Sigma = S))
}

(S <- covLomax(3, rep(1, 5))$Sigma)
(mu <- covLomax(3, rep(1, 5))$mu)

x <- NonNorMvtDist::rmvlomax(10000, parm1 = 3, parm2 = rep(1, 5))
cov(x)
colMeans(x)
