## --Filename: CF_fiveCopulas.R
######################### Miss rate ###########################################
library(copula)
source("Arc_copulas.R")
nx <- 500
p <- 5
five_Clayton <- claytonCopula(param = 2, dim = p)
X1 <- rCopula(500, five_Clayton)
five_Frank <- frankCopula(param = 2, dim = p)
five_Gumbel <-  gumbelCopula(param = 2, dim = p)
five_joe <- joeCopula(param = 2, dim = p)
J <- c(-1, 0, 1, 2)
set.seed(11)
A <- matrix(sample(J, size = p*p, replace = T), nrow = p, byrow = T)
diag(A)[diag(A) <= 0] <- rep(10, sum(diag(A)==0))
S <- t(A) %*% A
S <- diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
S <- round(S, digits = 1) # check with equal covariance. 
eigen(S,only.values = T) # Check S is pd 
five_t <- tCopula(S[upper.tri(S, diag = F)], dim = p, dispstr = "un",
                  df = 2, df.fixed = TRUE)

five_Normal <- normalCopula(S[upper.tri(S, diag = F)], dim = p, dispstr = "un")

myfiveCopulas <- list("Frank" = five_Frank, "Gumbel" = five_Gumbel, 
                      "Clayton" = five_Clayton, "Joe" = five_joe, 
                      "Normal" = five_Normal , "t" = five_t
)

#### Get misrate of Copulas classification using eigenstructure
myfiveRate <- list()
length(myfiveRate) <- length(myfiveCopulas)
names(myfiveRate) <- names(myfiveCopulas)
######################### Miss rate ###########################################
m <- 2e3; nx <- 1e3; myweight <- c(2/5, 1/5, 2/5) #change the number
print(m) 
print(p)
seed <- 1234; 
set.seed(seed)
for (i in 1:length(myfiveCopulas)){
  cops <- myfiveCopulas[[i]]
  print(i)
  re <- replicate(m, expr = {
    U <- copula::rCopula(nx, copula = cops)
    #### Get the estimate and generate reference samples
    frank.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::frankCopula(dim = p), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::frankCopula(dim = p), data = U, method = "itau"
        )
      })
    )
    frank.cop <- copula::frankCopula(param = frank.param, dim = p)
    frank.Y <- copula::rCopula(nx, copula = frank.cop)
    frank.Y2 <- copula::rCopula(2*nx, copula = frank.cop)
    frank.Y4 <- copula::rCopula(4*nx, copula = frank.cop)
    frank.Y8 <- copula::rCopula(8*nx, copula = frank.cop)
    frank.Y16 <- copula::rCopula(16*nx, copula = frank.cop)
    frank.Y32 <- copula::rCopula(32*nx, copula = frank.cop)
    
    gumbel.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::gumbelCopula(dim = p), data = U, method = "itau"
        )
      })
    )
    gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
    gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
    gumbel.Y2 <- copula::rCopula(2*nx, copula = gumbel.cop)
    gumbel.Y4 <- copula::rCopula(4*nx, copula = gumbel.cop)
    gumbel.Y8 <- copula::rCopula(8*nx, copula = gumbel.cop)
    gumbel.Y16 <- copula::rCopula(16*nx, copula = gumbel.cop)
    gumbel.Y32 <- copula::rCopula(32*nx, copula = gumbel.cop)
    
    clayton.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::claytonCopula(dim = p), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::claytonCopula(dim = p), data = U, method = "itau"
        )
      })
    )
    if (clayton.param < 0) {clayton.param <- 0}
    clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
    clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
    clayton.Y2 <- copula::rCopula(2*nx, copula = clayton.cop)
    clayton.Y4 <- copula::rCopula(4*nx, copula = clayton.cop)
    clayton.Y8 <- copula::rCopula(8*nx, copula = clayton.cop)
    clayton.Y16 <- copula::rCopula(16*nx, copula = clayton.cop)
    clayton.Y32 <- copula::rCopula(32*nx, copula = clayton.cop)
    
    joe.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::joeCopula(dim = p), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::joeCopula(dim = p), data = U, method = "itau"
        )
      })
    )
    joe.cop <- copula::joeCopula(param = joe.param, dim = p)
    joe.Y <- copula::rCopula(nx, copula = joe.cop)
    joe.Y2 <- copula::rCopula(2*nx, copula = joe.cop)
    joe.Y4 <- copula::rCopula(4*nx, copula = joe.cop)
    joe.Y8 <- copula::rCopula(8*nx, copula = joe.cop)
    joe.Y16 <- copula::rCopula(16*nx, copula = joe.cop)
    joe.Y32 <- copula::rCopula(32*nx, copula = joe.cop)
    
    normal.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::normalCopula(dim = p, dispstr = "un"),
        data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::normalCopula(dim = p, dispstr = "un"), 
          data = U, method = "itau"
        )
      })
    )
    normal.cop <- copula::normalCopula(param = normal.param, 
                                       dim = p, dispstr = "un"
    )
    normal.Y <- copula::rCopula(nx, copula = normal.cop)
    normal.Y2 <- copula::rCopula(2*nx, copula = normal.cop)
    normal.Y4 <- copula::rCopula(4*nx, copula = normal.cop)
    normal.Y8 <- copula::rCopula(8*nx, copula = normal.cop)
    normal.Y16 <- copula::rCopula(16*nx, copula = normal.cop)
    normal.Y32 <- copula::rCopula(32*nx, copula = normal.cop)
    
    t.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::tCopula(dim = p, df.fixed = FALSE, dispstr = "un"), 
        data = U, method = "mpl"
      ),
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::tCopula(dim = p, dispstr = "un"), 
          data = U, method = "itau.mpl"
        )
      })
    )
    t.cop <- copula::tCopula(param = t.param[-length(t.param)], dim = p, 
                             df = round(t.param["df"]), dispstr = "un"
    )
    t.Y <- copula::rCopula(nx, copula = t.cop)
    t.Y2 <- copula::rCopula(2*nx, copula = t.cop)
    t.Y4 <- copula::rCopula(4*nx, copula = t.cop)
    t.Y8 <- copula::rCopula(8*nx, copula = t.cop)
    t.Y16 <- copula::rCopula(16*nx, copula = t.cop)
    t.Y32 <- copula::rCopula(32*nx, copula = t.cop)
    
    lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y, normal.Y, t.Y)
    lstX2 <- list(frank.Y2, gumbel.Y2, clayton.Y2, joe.Y2, normal.Y2, t.Y2)
    lstX4 <- list(frank.Y4, gumbel.Y4, clayton.Y4, joe.Y4, normal.Y4, t.Y4)
    lstX8 <- list(frank.Y8, gumbel.Y8, clayton.Y8, joe.Y8, normal.Y8, t.Y8)
    lstX16 <- list(frank.Y16, gumbel.Y16, clayton.Y16, 
                   joe.Y16, normal.Y16, t.Y16
    )
    lstX32 <- list(frank.Y32, gumbel.Y32, clayton.Y32, 
                   joe.Y32, normal.Y32, t.Y32
    )
    lstCop <- list(frank.cop, gumbel.cop, clayton.cop, 
                   joe.cop, normal.cop, t.cop
    )
    
    c(
      Eigen = mul.discr.eigen(lstX, U)$pop, 
      Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop,
      Eigen_reg2 = mul.discr.eigen.reg(lstX2, U, w = rep(1, 3))$pop,
      Eigen_reg4 = mul.discr.eigen.reg(lstX4, U, w = rep(1, 3))$pop,
      Eigen_reg8 = mul.discr.eigen.reg(lstX8, U, w = rep(1, 3))$pop,
      Eigen_reg16 = mul.discr.eigen.reg(lstX16, U, w = rep(1, 3))$pop,
      Eigen_reg32 = mul.discr.eigen.reg(lstX32, U, w = rep(1, 3))$pop,
      Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
      NP = mul.NP(lstCop, U)$pop)
  })
  
  five_rate_tmp <- t(
    apply(re, MARGIN = 1, 
          FUN =  function(arr) {
            c(mean(arr== "Population 1"), 
              mean(arr=="Population 2"),
              mean(arr == "Population 3"),
              mean(arr== "Population 4"), 
              mean(arr== "Population 5"), 
              mean(arr== "Population 6")
            )
          }
    )
  )
  colnames(five_rate_tmp) <- c("Frank", "Gumbel", "Clayton", "Joe", "Normal", "t")
  # colnames(five_rate_tmp) <- c("Frank", "Gumbel", "Clayton", "Joe", "Normal")
  myfiveRate[[i]] <- five_rate_tmp
}
myfiveRate
tmp <- myfiveRate[[1]][, 1]
for (k in 2:6){
  tmp <- cbind(tmp, myfiveRate[[k]][, k])
}
tmp
xtable::xtable(tmp, digits = 4)
tmp <- myfiveRate[[1]]

for (k in 1:9){
  tmp <- sapply(1:6, function(x) myfiveRate[[x]][k, ])
  diag(tmp) <- NA
  tmp1 <- rbind(apply(tmp, MARGIN = 2, function(x) which.max(x)),
                apply(tmp, MARGIN = 2, function(x) max(x, na.rm = T))
  )
  colnames(tmp1) <- c("Frank", "Gumbel", "Clayton", "Joe", "Normal", "t")
  tmp1[1, ] <- dplyr::case_when(tmp1[1, ] == 1 ~ "Frank",
                                tmp1[1, ] == 2 ~ "Gumbel",
                                tmp1[1, ] == 3 ~ "Clayton",
                                tmp1[1, ] == 4 ~ "Joe",
                                tmp1[1, ] == 5 ~ "Normal",
                                tmp1[1, ] == 6 ~ "t",
  )
  
  tmp1 <- xtable::xtable(tmp1, digits = 4)
  print(tmp1)
}


