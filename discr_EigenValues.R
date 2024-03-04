## -- Filename: discr_Eigenvalues.R
## Group discrimination using eigenvalue 
antieigen <- function(X, U, type = 1){
  # Methods calculating eigenvalues
  # eta1 - H matrix
  p <- ncol(U)
  n <- nrow(U)
  r <- nrow(X)
  if (type  == 1){
    H <- (n + r)/r * solve((t(X)%*% X + t(U) %*% U )) %*% (t(X) %*% X)
    mye <- Re(eigen(H)$values)
    mu <- (2 * sqrt(mye[1] * mye[p])) / (mye[1] + mye[p])
  }
  
  if (type == 2){
    H <- (n + r)/r * solve((t(X)%*% X + t(U) %*% U )) %*% (t(X) %*% X)
    mye <- Re(eigen(H)$values)
    mu <- prod(
      sapply(
        1:ceiling(p/2),
        function(i) (2 * sqrt(mye[i] * mye[p-i+1]))/(mye[i] + mye[p - i +1])
        )
      )
  }

  if (type == 3){
    D <- n/r * solve(t(U) %*% U) %*% (t(X) %*% X)
    mye <- Re(eigen(D)$values)
    mu <- (2 * sqrt(mye[1] * mye[p])) / (mye[1] + mye[p])
  }
  
  if (type == 4){
    D <- n/r * solve(t(U) %*% U) %*% (t(X) %*% X)
    mye <- Re(eigen(D)$values)
    mu <- prod(
      sapply(
        1:ceiling(p/2),
        function(i) (2 * sqrt(mye[i] * mye[p-i+1]))/(mye[i] + mye[p - i +1])
        )
      )
  }
  return(mu)
}
#####################
discr.antieigen <- function(X1, X2, U, type = 1){
  mu1 <- antieigen(X1, U, type = type)
  mu2 <- antieigen(X2, U, type = type)
  deci <- ifelse(mu1 == mu2, "Can not decide", 
                 ifelse(mu1 > mu2, "Population 1", "Population 2")
                 )
  return(list(pop_1 = mu1, pop_2 = mu2, pop = deci))
}
############################################
############################################
antieigen.reg <- function(X, U, cut.point = c(.3, .7)){
  #calculate antieigen values bX each region
  re <- c(
    antieigen(X[which(X[, 1] <=cut.point[1]), ], 
              U[which(U[,1] <= cut.point[1]), ]
    ), 
    antieigen(X[which((X[, 1] > cut.point[1])&(X[, 1] <=cut.point[2])), ],
              U[which((U[, 1] > cut.point[1])&(U[, 1] <=cut.point[2])), ]
              ), 
    antieigen(X[which(X[, 1] > cut.point[2]), ],
              U[which(U[,1] > cut.point[2]), ]
              )
    )
  return(re) #the higher the more similar
}
