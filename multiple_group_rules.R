## --Filename: multiple_group_rules
# lstX: list of references samples
# U: "to be classified" sample
#############################
mul.discr.linear <- function(lstX, U){
  require(MASS)
  k <- length(lstX)
  DF <- lapply(1:k, FUN = function(j) {
    df <- data.frame(lstX[[j]])
    df$pop <- j
    df}
    )
  DF <- Reduce(rbind, DF)
  DF$pop <- as.factor(DF$pop)
  lda_model <- MASS::lda(pop ~., data = DF)
  
  lda_predictions <- predict(lda_model, data.frame(U))
  re <- table(lda_predictions$class)
  deci <- paste("Population", as.character(which.max(re)))
  # deci <- ifelse(num1 >= nU/2, "Population 1", "Population 2" )
  return(list(pop = deci, value = re))
}
###############################
###############################
mul.discr.MaxSep <- function(lstX, U){
  require(MASS)
  p <- ncol(U); nU <- nrow(U)
  k <- length(lstX)
  DF <- lapply(1:k, FUN = function(j) {
    df <- data.frame(lstX[[j]])
    df$pop <- j
    df}
  )
  DF <- Reduce(rbind, DF)
  DF$pop <- as.factor(DF$pop)
  lda_model <- MASS::lda(pop ~., data = DF)
  
  bigD <- array(NA, dim = c(nrow(U), length(lstX))) # For U to each pop
  lda_U <- U %*% lda_model$scaling
  for (i in 1:length(lstX)){
    lda_functions <- as.matrix(
      subset(DF, pop == as.character(i), select = -(p+1)) 
    ) %*% lda_model$scaling
    d <- apply(
      lda_functions, MARGIN = 1, 
      FUN = function(u) mean(norm(lda_U - rep(u, nU), "F"))
    )
    bigD[, i] <- d
  }
  
  re <- table(apply(bigD, MARGIN = 1, which.min))
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = re))
}
###############################
###############################
mul.discr.averageD <- function(lstX, U){
  nU <- nrow(U)
  bigD <- array(NA, dim = c(nrow(U), length(lstX)))
  for (i in 1:length(lstX)){
    X <- lstX[[i]]
    invX <- solve(cov(X))
    nX <- nrow(X)
    d <- apply(U, MARGIN = 1, FUN = function(x){
      mean(
        sqrt(
          unlist(
            apply(
              X, MARGIN = 1, FUN = function(y) sqrt((x - y) %*% invX %*% (x - y))
              )
            )
          )
        )
    })
    bigD[, i] <- d
  }
  re <- table(apply(bigD, MARGIN = 1, which.min))
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = re))
}
###############################
###############################
mul.discr.leverage <- function(lstX, U){
  nU <- nrow(U)
  bigH <- array(NA, dim = c(nrow(U), length(lstX)))
  for  (i in 1:length(lstX)){
    X <- lstX[[i]]
    invX <- solve(t(X) %*% X)
    h <- apply(U, MARGIN = 1, FUN = function(x) t(x) %*% invX %*% x )
    bigH[, i] <- h
  }
  re <- table(apply(bigH, MARGIN = 1, which.min))
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = re))
}
###############################
###############################
mul.discr.eigen <- function(lstX, U, type = 1){
  re <- lapply(lstX,antieigen, U = U, type = type)  
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = unlist(re)))
}
###############################
###############################
mul.discr.eigen.reg <- function(lstX, U, 
                                cut.point = c(.35, .65), 
                                w = rep(1, length(cut.point) + 1)
                                ){
 # w is the weight 
  re <- lapply(
    lstX, FUN = function(x) weighted.mean(
      antieigen.reg(x, U,  cut.point = cut.point), w = w
    )
  )
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = unlist(re)))
}
###############################
###############################
mul.Emp.copulas <- function(lstCop, U){
  require(copula)
  ############################
  eC <- emp.copulas(U)
  nU <- nrow(U)
  ############################
  myfunc <- function(Ck){
    eCk <- matrix(nrow = nx, ncol = nx)
    for (i in 1:nU){
      for (j in 1:nU){
        eCk[i,j] <- copula::pCopula(u = c(i/nx, j/nx), copula = Ck)
      }
    }
    re <- norm(eC- eCk, "F")
    return(re)
  }
  ############################
  re <- c()
  for (k in 1:length(lstCop)){
    Ck <- lstCop[[k]]
    re <- c(re, myfunc(Ck))
  }
  
  deci <- paste("Population", as.character(which.min(re)))
  return(list(pop = deci, value = re))
}
###############################
###############################
mul.NP <- function(lstCop, U){
  require(copula)
  re <- c()
  ###############################
  for (k in 1:length(lstCop)){
    Ck <- lstCop[[k]]
    re_tmp <- sum(
      apply(
        U, MARGIN = 1, function(u) copula::dCopula(u, copula = Ck, log = T)
        )
      )
    re <- c(re, re_tmp)
  }
  deci <- paste("Population", as.character(which.max(re)))
  return(list(pop = deci, value = re))
}
