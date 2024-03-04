emp.copulas <- function(X){
  #Emperical copulas is not unique
  #calculate the Deheuvels emperical bivariate copulas
  nx <- nrow(X)
  
  sortX <- apply(X, MARGIN = 2, FUN = sort)
  eC <- matrix(0,nrow = nx, ncol= nx)
  
  for (i in 1:nx){
    for (j in 1:nx){ 
      for (k in 1:nx){
        if (X[k, 1] <= sortX[i, 1] & X[k, 2] <= sortX[j, 2] )
          eC[i, j] <- eC[i, j] + 1
      }
    }
  }
  return(eC/nx)
}

# X <- matrix(c(0.95, 0.24,
#               0.53, 0.16, 
#               0.77, 0.56,
#               0.19, 0.33,
#               0.32, 0.8), ncol = 2, byrow = T)
# emp.copulas(X)

######################
######################
Ecopulas.distance <- function(X, Ck){
  require(copula)
  # calculate the distance of Eperical copulas and hypothesized copulas 
  eC <- emp.copulas(X)
  nx <- nrow(X)
  eCk <- matrix(nrow = nx, ncol = nx)
  for (i in 1:nx){
    for (j in 1:nx){
      eCk[i,j] <- copula::pCopula(u = c(i/nx, j/nx), copula = Ck)
    }
  }
  re <- norm((eC - eCk), "F")
  # the smalleer the better
  return(re)
}

# Ecopulas.distance(X.frank, Ck = frank.cop)
