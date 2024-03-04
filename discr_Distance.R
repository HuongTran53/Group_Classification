## --Filename: discr_Distance.R
## Group discriminations my majority classification rules 
########################################
discr.linear <- function(X1 = NULL, X2 = NULL, U, true.params = T, 
                         params = ifelse(
                           true.params, list(mu1, mu2, Sigma), NULL)
                         ){
  nU <- nrow(U)
  if (!true.params){
    warning("Unknown true parameter, estimators will be uses")
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    mu1 <- colMeans(X1)
    mu2 <- colMeans(X2)
    
    S1 <- cov(X1)
    S2 <- cov(X2)
    Sigma <- ((n1 - 1) * S1 + (n2 - 1) *S2)/(n1 + n2 - 2)
  } else {
    warning("Using true parameters")
    mu1 <- params$mu1
    mu2 <- params$mu2
    Sigma <- params$Sigma
  }
  # assume that the two training set have equal variance
  a <- solve(Sigma) %*% (mu1 - mu2)
  myu <- U %*% a - as.numeric(.5 * t(a) %*% (mu1 + mu2))
  re <- ifelse(myu >= 0, "Population 1", "Population 2")
  num1 <- sum(re == "Population 1")
  deci <- ifelse(num1 == nU/2, "Can not decide", 
                 ifelse(num1 > nU/2, "Population 1", "Population 2" ))
  return(list(pop_1 = num1, pop_2 = nU - num1, pop = deci))
}
########################################
########################################
discr.MaxSep <- function(X1 = NULL, X2 = NULL, U, true.params = T, 
                         params = ifelse(
                           true.params, list(mu1, mu2, Sigma), NULL)
                         ){
  # apply linear discriminant
  # average them by each available sample before classification 
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  nU <- nrow(U)
  if (true.params){
    warning("Using true parameters")
    mu1 <- params$mu1
    mu2 <- params$mu2
    Sigma <- params$Sigma
    a <- solve(Sigma) %*% (mu1 - mu2) 
    bigU <- apply(U,1 , 
                  function(y) t(a) %*% ((y - mu1) %*% t((y - mu1)) - 
                                          (y - mu2)%*% t((y - mu2)))%*% a
                  )
    re <- ifelse(bigU <= 0, "Population 1", "Population 2")
  } else { 
    warning("Unknown true parameter, estimators will be uses")
    mu1 <- colMeans(X1)
    mu2 <- colMeans(X2)
    S1 <- cov(X1)
    S2 <- cov(X2)
    Sigma <- ((n1 - 1) * S1 + (n2 - 1) *S2)/(n1 + n2 - 2)
    a <- solve(Sigma) %*% (mu1 - mu2) 
    w1 <- X1 %*% a
    w2 <- X2 %*% a
    dU <- U %*% a
    bigU <- apply(dU, 1, function(x) mean((w1 - x))^2 - mean((w2 - x))^2)
    re <- ifelse(bigU <= 0, "Population 1", "Population 2")
  }
  num1 <- sum(re == "Population 1")
  deci <- ifelse(num1 == nU/2, "Can not decide", 
                 ifelse(num1 > nU/2, "Population 1", "Population 2" )
                 )
  return(list(pop_1 = num1, pop_2 = nU - num1, pop = deci))
}
########################################
########################################
discr.averageD <- function(X1, X2,  U){
  invX1 <- solve(cov(X1))
  invX2 <- solve(cov(X2))
  n1 <- nrow(X1); n2 <- nrow(X2); nU <- nrow(U)
  d1 <- apply(U, MARGIN = 1, FUN = function(x){
    mean(sqrt(unlist(apply(X1, MARGIN = 1, 
                           FUN = function(y) sqrt((x - y) %*% invX1 %*% (x - y))
                           )
                     )
              )
         )
  }
  )
  d2 <- apply(U, MARGIN = 1, FUN = function(x){
    mean(sqrt(unlist(apply(X2, MARGIN = 1, 
                           FUN = function(y) sqrt((x - y) %*% invX2 %*% (x - y))
                           )
                     )
              )
         )
  }
  )
  re <- ifelse(d1 <= d2, "Population 1", "Population 2")
  num1 <- sum(re == "Population 1")
  deci <- ifelse(num1 == nU/2, "Can not decide", 
                 ifelse(num1 > nU/2, "Population 1", "Population 2" )
                 )
  return(list(pop_1 = num1, pop_2 = nU - num1, pop = deci))
}
########################################
########################################
discr.leverage <- function(X1, X2, U){
  nU <- nrow(U)
  n1 <- nrow(X1)
  n2 <- nrow(X1)
  invX1 <- (n1 - 1)*solve(t(X1) %*% X1)
  invX2 <- (n2 - 1)*solve(t(X2) %*% X2)
  H1 <- apply(U, MARGIN = 1, FUN = function(x) t(x) %*% invX1 %*% x )
  H2 <- apply(U, MARGIN = 1, FUN = function(x) t(x) %*% invX2 %*% x)
  re <- ifelse(H1 <= H2, "Population 1", "Population 2")
  num1 <- sum(re == "Population 1")
  deci <- ifelse(num1 == nU/2, "Can not decide", 
                 ifelse(num1 > nU/2, "Population 1", "Population 2" )
                 )
  return(list(pop_1 = num1, pop_2 = nU - num1, pop = deci))
}
