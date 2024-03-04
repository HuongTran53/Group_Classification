# # Get copulas data generator 
require(copula)
source("multiple_group_rules.R")
source("discr_Distance.R")
source("discr_EigenValues.R")
source("emp_copulas.R")
source("Arc_copulas.R")


######################### Miss rate ###########################################
#save(biFrank_rate, file = "biFrank_rate.RDATA")
# Gumbel
set.seed(1234)
re2 <- replicate(m, expr = {

  U <- copula::rCopula(nx, copula = bi_Gumbel)
  
  frank.param <- coef(copula::fitCopula(copula = copula::frankCopula(dim = p), data = U, method = "mpl"))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  frank.Y <- copula::rCopula(nx, copula = frank.cop)
  frank.Y2 <- copula::rCopula(2*nx, copula = frank.cop)
  frank.Y4 <- copula::rCopula(4*nx, copula = frank.cop)
  frank.Y8 <- copula::rCopula(8*nx, copula = frank.cop)
  frank.Y16 <- copula::rCopula(16*nx, copula = frank.cop)
  frank.Y32 <- copula::rCopula(32*nx, copula = frank.cop)
  # frank.Y100 <- copula::rCopula(100*nx, copula = frank.cop)
  
  gumbel.param <-coef(copula::fitCopula(copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
  gumbel.Y2 <- copula::rCopula(2*nx, copula = gumbel.cop)
  gumbel.Y4 <- copula::rCopula(4*nx, copula = gumbel.cop)
  gumbel.Y8 <- copula::rCopula(8*nx, copula = gumbel.cop)
  gumbel.Y16 <- copula::rCopula(16*nx, copula = gumbel.cop)
  gumbel.Y32 <- copula::rCopula(32*nx, copula = gumbel.cop)
  # gumbel.Y100 <- copula::rCopula(100*nx, copula = gumbel.cop)
  
  clayton.param <- coef(copula::fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "mpl"))
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
  clayton.Y2 <- copula::rCopula(2*nx, copula = clayton.cop)
  clayton.Y4 <- copula::rCopula(4*nx, copula = clayton.cop)
  clayton.Y8 <- copula::rCopula(8*nx, copula = clayton.cop)
  clayton.Y16 <- copula::rCopula(16*nx, copula = clayton.cop)
  clayton.Y32 <- copula::rCopula(32*nx, copula = clayton.cop)
  # clayton.Y100 <- copula::rCopula(100*nx, copula = clayton.cop)
  
  joe.param <- coef(copula::fitCopula(copula = copula::joeCopula(dim= p), data = U, method = "mpl"))
  joe.cop <- copula::joeCopula(param = joe.param, dim = p)
  joe.Y <-  copula::rCopula(nx, copula = joe.cop)
  joe.Y2 <-  copula::rCopula(2*nx, copula = joe.cop)
  joe.Y4 <-  copula::rCopula(4*nx, copula = joe.cop)
  joe.Y8 <-  copula::rCopula(8*nx, copula = joe.cop)
  joe.Y16 <-  copula::rCopula(16*nx, copula = joe.cop)
  joe.Y32 <-  copula::rCopula(32*nx, copula = joe.cop)
  # joe.Y100 <-  copula::rCopula(100*nx*p, copula = joe.cop)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y)
  lstX2 <- list(frank.Y2, gumbel.Y2, clayton.Y2, joe.Y2)
  lstX4 <- list(frank.Y4, gumbel.Y4, clayton.Y4, joe.Y4)
  lstX8 <- list(frank.Y8, gumbel.Y8, clayton.Y8, joe.Y8)
  lstX16 <- list(frank.Y16, gumbel.Y16, clayton.Y16, joe.Y16)
  lstX32 <- list(frank.Y32, gumbel.Y32, clayton.Y32, joe.Y32)
  # lstX100 <- list(frank.Y100, gumbel.Y100, clayton.Y100, joe.Y100)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, joe.cop)
  
  
  c(
    discr_linear = mul.discr.linear(lstX, U)$pop,
    MaxSep = mul.discr.MaxSep(lstX, U)$pop,
    AverageD = mul.discr.averageD(lstX, U)$pop,
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop,
    Eigen_reg2 = mul.discr.eigen.reg(lstX2, U, w = rep(1, 3))$pop,
    Eigen_reg4 = mul.discr.eigen.reg(lstX4, U, w = rep(1, 3))$pop,
    Eigen_reg8 = mul.discr.eigen.reg(lstX8, U, w = rep(1, 3))$pop,
    Eigen_reg16 = mul.discr.eigen.reg(lstX16, U, w = rep(1, 3))$pop,
    Eigen_reg32 = mul.discr.eigen.reg(lstX32, U, w = rep(1, 3))$pop,
    # Eigen_reg100 = mul.discr.eigen.reg(lstX100, U, w = rep(1, 3))$pop,
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)

})
biGumbel_rate <- t(apply(re2, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(biGumbel_rate) <- c("Frank", "Gumbel", "Clayton", "Joe")
#save(biGumbel_rate, file = "biGumbel_rate.RDATA")
# Clayton
set.seed(1234)
re3 <- replicate(m, expr = {
  U <- copula::rCopula(nx, copula = bi_Clayton)

  frank.param <- coef(copula::fitCopula(copula = copula::frankCopula(dim = p), data = U, method = "mpl"))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  frank.Y <- copula::rCopula(nx, copula = frank.cop)
  frank.Y2 <- copula::rCopula(2*nx, copula = frank.cop)
  frank.Y4 <- copula::rCopula(4*nx, copula = frank.cop)
  frank.Y8 <- copula::rCopula(8*nx, copula = frank.cop)
  frank.Y16 <- copula::rCopula(16*nx, copula = frank.cop)
  frank.Y32 <- copula::rCopula(32*nx, copula = frank.cop)
  # frank.Y100 <- copula::rCopula(100*nx, copula = frank.cop)
  
  gumbel.param <-coef(copula::fitCopula(copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
  gumbel.Y2 <- copula::rCopula(2*nx, copula = gumbel.cop)
  gumbel.Y4 <- copula::rCopula(4*nx, copula = gumbel.cop)
  gumbel.Y8 <- copula::rCopula(8*nx, copula = gumbel.cop)
  gumbel.Y16 <- copula::rCopula(16*nx, copula = gumbel.cop)
  gumbel.Y32 <- copula::rCopula(32*nx, copula = gumbel.cop)
  # gumbel.Y100 <- copula::rCopula(100*nx, copula = gumbel.cop)
  
  clayton.param <- coef(tryCatch(
    fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "mpl"), 
    error = function(ex){
      warning("MLE method has error - ", ex)
      copula::fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "itau")
    }))
  
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
  clayton.Y2 <- copula::rCopula(2*nx, copula = clayton.cop)
  clayton.Y4 <- copula::rCopula(4*nx, copula = clayton.cop)
  clayton.Y8 <- copula::rCopula(8*nx, copula = clayton.cop)
  clayton.Y16 <- copula::rCopula(16*nx, copula = clayton.cop)
  clayton.Y32 <- copula::rCopula(32*nx, copula = clayton.cop)
  # clayton.Y100 <- copula::rCopula(100*nx, copula = clayton.cop)
  
  joe.param <- coef(copula::fitCopula(copula = copula::joeCopula(dim= p), data = U, method = "mpl"))
  joe.cop <- copula::joeCopula(param = joe.param, dim = p)
  joe.Y <-  copula::rCopula(nx, copula = joe.cop)
  joe.Y2 <-  copula::rCopula(2*nx, copula = joe.cop)
  joe.Y4 <-  copula::rCopula(4*nx, copula = joe.cop)
  joe.Y8 <-  copula::rCopula(8*nx, copula = joe.cop)
  joe.Y16 <-  copula::rCopula(16*nx, copula = joe.cop)
  joe.Y32 <-  copula::rCopula(32*nx, copula = joe.cop)
  # joe.Y100 <-  copula::rCopula(100*nx*p, copula = joe.cop)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y)
  lstX2 <- list(frank.Y2, gumbel.Y2, clayton.Y2, joe.Y2)
  lstX4 <- list(frank.Y4, gumbel.Y4, clayton.Y4, joe.Y4)
  lstX8 <- list(frank.Y8, gumbel.Y8, clayton.Y8, joe.Y8)
  lstX16 <- list(frank.Y16, gumbel.Y16, clayton.Y16, joe.Y16)
  lstX32 <- list(frank.Y32, gumbel.Y32, clayton.Y32, joe.Y32)
  # lstX100 <- list(frank.Y100, gumbel.Y100, clayton.Y100, joe.Y100)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, joe.cop)
  
  
  c(
    discr_linear = mul.discr.linear(lstX, U)$pop,
    MaxSep = mul.discr.MaxSep(lstX, U)$pop,
    AverageD = mul.discr.averageD(lstX, U)$pop,
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop,
    Eigen_reg2 = mul.discr.eigen.reg(lstX2, U, w = rep(1, 3))$pop,
    Eigen_reg4 = mul.discr.eigen.reg(lstX4, U, w = rep(1, 3))$pop,
    Eigen_reg8 = mul.discr.eigen.reg(lstX8, U, w = rep(1, 3))$pop,
    Eigen_reg16 = mul.discr.eigen.reg(lstX16, U, w = rep(1, 3))$pop,
    Eigen_reg32 = mul.discr.eigen.reg(lstX32, U, w = rep(1, 3))$pop,
    # Eigen_reg100 = mul.discr.eigen.reg(lstX100, U, w = rep(1, 3))$pop,
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)
  
})
biClayton_rate <- t(apply(re3, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(biClayton_rate) <- c("Frank", "Gumbel", "Clayton", "Joe")
#save(biClayton_rate, file = "biClayton_rate.RDATA")
# joe
set.seed(1234)
re4 <- replicate(m, expr = {

  U <- copula::rCopula(nx, copula = bi_joe)

  frank.param <- coef(copula::fitCopula(copula = copula::frankCopula(dim = p), data = U, method = "mpl"))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  frank.Y <- copula::rCopula(nx, copula = frank.cop)
  frank.Y2 <- copula::rCopula(2*nx, copula = frank.cop)
  frank.Y4 <- copula::rCopula(4*nx, copula = frank.cop)
  frank.Y8 <- copula::rCopula(8*nx, copula = frank.cop)
  frank.Y16 <- copula::rCopula(16*nx, copula = frank.cop)
  frank.Y32 <- copula::rCopula(32*nx, copula = frank.cop)
  # frank.Y100 <- copula::rCopula(100*nx, copula = frank.cop)
  
  gumbel.param <-coef(copula::fitCopula(copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
  gumbel.Y2 <- copula::rCopula(2*nx, copula = gumbel.cop)
  gumbel.Y4 <- copula::rCopula(4*nx, copula = gumbel.cop)
  gumbel.Y8 <- copula::rCopula(8*nx, copula = gumbel.cop)
  gumbel.Y16 <- copula::rCopula(16*nx, copula = gumbel.cop)
  gumbel.Y32 <- copula::rCopula(32*nx, copula = gumbel.cop)
  gumbel.Y100 <- copula::rCopula(100*nx, copula = gumbel.cop)
  
  clayton.param <- coef(tryCatch(
    fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "mpl"), 
    error = function(ex){
      warning("MLE method has error - ", ex)
      copula::fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "itau")
    }))
  
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
  clayton.Y2 <- copula::rCopula(2*nx, copula = clayton.cop)
  clayton.Y4 <- copula::rCopula(4*nx, copula = clayton.cop)
  clayton.Y8 <- copula::rCopula(8*nx, copula = clayton.cop)
  clayton.Y16 <- copula::rCopula(16*nx, copula = clayton.cop)
  clayton.Y32 <- copula::rCopula(32*nx, copula = clayton.cop)
  # clayton.Y100 <- copula::rCopula(100*nx, copula = clayton.cop)
  
  joe.param <- coef(copula::fitCopula(copula = copula::joeCopula(dim= p), data = U, method = "mpl"))
  joe.cop <- copula::joeCopula(param = joe.param, dim = p)
  joe.Y <-  copula::rCopula(nx, copula = joe.cop)
  joe.Y2 <-  copula::rCopula(2*nx, copula = joe.cop)
  joe.Y4 <-  copula::rCopula(4*nx, copula = joe.cop)
  joe.Y8 <-  copula::rCopula(8*nx, copula = joe.cop)
  joe.Y16 <-  copula::rCopula(16*nx, copula = joe.cop)
  joe.Y32 <-  copula::rCopula(32*nx, copula = joe.cop)
  # joe.Y100 <-  copula::rCopula(100*nx*p, copula = joe.cop)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y)
  lstX2 <- list(frank.Y2, gumbel.Y2, clayton.Y2, joe.Y2)
  lstX4 <- list(frank.Y4, gumbel.Y4, clayton.Y4, joe.Y4)
  lstX8 <- list(frank.Y8, gumbel.Y8, clayton.Y8, joe.Y8)
  lstX16 <- list(frank.Y16, gumbel.Y16, clayton.Y16, joe.Y16)
  lstX32 <- list(frank.Y32, gumbel.Y32, clayton.Y32, joe.Y32)
  # lstX100 <- list(frank.Y100, gumbel.Y100, clayton.Y100, joe.Y100)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, joe.cop)
  
  
  c(
    discr_linear = mul.discr.linear(lstX, U)$pop,
    MaxSep = mul.discr.MaxSep(lstX, U)$pop,
    AverageD = mul.discr.averageD(lstX, U)$pop,
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop,
    Eigen_reg2 = mul.discr.eigen.reg(lstX2, U, w = rep(1, 3))$pop,
    Eigen_reg4 = mul.discr.eigen.reg(lstX4, U, w = rep(1, 3))$pop,
    Eigen_reg8 = mul.discr.eigen.reg(lstX8, U, w = rep(1, 3))$pop,
    Eigen_reg16 = mul.discr.eigen.reg(lstX16, U, w = rep(1, 3))$pop,
    Eigen_reg32 = mul.discr.eigen.reg(lstX32, U, w = rep(1, 3))$pop,
    # Eigen_reg100 = mul.discr.eigen.reg(lstX100, U, w = rep(1, 3))$pop,
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)

})

bijoe_rate <- t(apply(re4, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(bijoe_rate) <- c("Frank", "Gumbel", "Clayton", "Joe")
#save(bijoe_rate, file = "bijoe_rate.RDATA")
##########################################
##########################################
# xtable::xtable(data.frame(biClayton_rate))
biFrank_rate
biClayton_rate
biGumbel_rate
bijoe_rate
