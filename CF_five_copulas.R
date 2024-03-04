######################### Miss rate ###########################################
tau <- 0.5; p <- 5
five_Clayton <- Arc_copula(p = p, tau = tau)$clayton.cop
five_Frank <- Arc_copula(p = p, tau = tau)$frank.cop
five_Gumbel <-  Arc_copula(p = p, tau = tau)$gumbel.cop
five_Normal <- Arc_copula(p = p, tau = tau)$normal.cop

######################### Miss rate ###########################################
m <- 1e3; nx <- 200; myweight <- c(2/5, 1/5, 2/5) #change the number
seed <- 1234; ny <- 2*nx
# Frank
set.seed(seed)
re1 <- replicate(m, expr = {
  
  U <- copula::rCopula(nx, copula = five_Frank)
  
  frank.param <- coef(copula::fitCopula(copula = copula::frankCopula(dim = p), data = U, method = "mpl"))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  frank.Y <- copula::rCopula(nx, copula = frank.cop)
  frank.Y2 <- copula::rCopula(2*nx, copula = frank.cop)
  frank.Y4 <- copula::rCopula(4*nx, copula = frank.cop)
  frank.Y8 <- copula::rCopula(8*nx, copula = frank.cop)
  frank.Y16 <- copula::rCopula(16*nx, copula = frank.cop)
  frank.Y32 <- copula::rCopula(32*nx, copula = frank.cop)
  frank.Y100 <- copula::rCopula(100*nx, copula = frank.cop)
  
  gumbel.param <-coef(copula::fitCopula(copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
  gumbel.Y2 <- copula::rCopula(2*nx, copula = gumbel.cop)
  gumbel.Y4 <- copula::rCopula(4*nx, copula = gumbel.cop)
  gumbel.Y8 <- copula::rCopula(8*nx, copula = gumbel.cop)
  gumbel.Y16 <- copula::rCopula(16*nx, copula = gumbel.cop)
  gumbel.Y32 <- copula::rCopula(32*nx, copula = gumbel.cop)
  gumbel.Y100 <- copula::rCopula(100*nx, copula = gumbel.cop)
  
  clayton.param <- coef(copula::fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "mpl"))
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
  clayton.Y2 <- copula::rCopula(2*nx, copula = clayton.cop)
  clayton.Y4 <- copula::rCopula(4*nx, copula = clayton.cop)
  clayton.Y8 <- copula::rCopula(8*nx, copula = clayton.cop)
  clayton.Y16 <- copula::rCopula(16*nx, copula = clayton.cop)
  clayton.Y32 <- copula::rCopula(32*nx, copula = clayton.cop)
  clayton.Y100 <- copula::rCopula(100*nx, copula = clayton.cop)
  
  normal.param <- coef(copula::fitCopula(copula = copula::normalCopula(dim= p), data = U, method = "mpl"))
  normal.cop <- copula::normalCopula(param = normal.param, dim = p)
  normal.Y <-  copula::rCopula(nx, copula = normal.cop)
  normal.Y2 <-  copula::rCopula(2*nx, copula = normal.cop)
  normal.Y4 <-  copula::rCopula(4*nx, copula = normal.cop)
  normal.Y8 <-  copula::rCopula(8*nx, copula = normal.cop)
  normal.Y16 <-  copula::rCopula(16*nx, copula = normal.cop)
  normal.Y32 <-  copula::rCopula(32*nx, copula = normal.cop)
  normal.Y100 <-  copula::rCopula(100*nx*p, copula = normal.cop)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, normal.Y)
  lstX2 <- list(frank.Y2, gumbel.Y2, clayton.Y2, normal.Y2)
  lstX4 <- list(frank.Y4, gumbel.Y4, clayton.Y4, normal.Y4)
  lstX8 <- list(frank.Y8, gumbel.Y8, clayton.Y8, normal.Y8)
  lstX16 <- list(frank.Y16, gumbel.Y16, clayton.Y16, normal.Y16)
  lstX32 <- list(frank.Y32, gumbel.Y32, clayton.Y32, normal.Y32)
  lstX100 <- list(frank.Y100, gumbel.Y100, clayton.Y100, normal.Y100)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, normal.cop)
  
  
  c(
    # discr_linear = mul.discr.linear(lstX, U)$pop, 
    # MaxSep = mul.discr.MaxSep(lstX, U)$pop, 
    # AverageD = mul.discr.averageD(lstX, U)$pop, 
    # Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop,
    Eigen_reg2 = mul.discr.eigen.reg(lstX2, U, w = rep(1, 3))$pop,
    Eigen_reg4 = mul.discr.eigen.reg(lstX4, U, w = rep(1, 3))$pop,
    Eigen_reg8 = mul.discr.eigen.reg(lstX8, U, w = rep(1, 3))$pop,
    Eigen_reg16 = mul.discr.eigen.reg(lstX16, U, w = rep(1, 3))$pop,
    Eigen_reg32 = mul.discr.eigen.reg(lstX32, U, w = rep(1, 3))$pop,
    Eigen_reg100 = mul.discr.eigen.reg(lstX100, U, w = rep(1, 3))$pop,
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)
})

fiveFrank_rate <- t(apply(re1, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(fiveFrank_rate) <- c("Frank", "Gumbel", "Clayton", "Normal")

#########################
# Gumbel
set.seed(seed)
re2 <- replicate(m, expr = {
  
  U <- copula::rCopula(nx, copula = five_Gumbel)
  
  frank.param <- coef(copula::fitCopula(copula = copula::frankCopula(dim = p), data = U, method = "mpl"))
  frank.Y <- copula::rCopula(nx, copula = copula::frankCopula(param = frank.param, dim = p))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  
  gumbel.param <-coef(copula::fitCopula(copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.Y <- copula::rCopula(nx, copula = copula::gumbelCopula(param = gumbel.param, dim = p))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  
  
  clayton.param <- coef(copula::fitCopula(copula = copula::claytonCopula(dim = p), data = U, method = "mpl"))
  clayton.Y <- copula::rCopula(nx, copula = copula::claytonCopula(param = clayton.param, dim = p))
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  
  normal.param <- coef(copula::fitCopula(copula = normalCopula(dim= p), data = U, method = "mpl"))
  normal.Y <- copula::rCopula(nx, copula = normalCopula(param = normal.param, dim = p))
  normal.cop <- copula::normalCopula(param = normal.param, dim = p)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, normal.Y)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, normal.cop)
  
  c(discr_linear = mul.discr.linear(lstX, U)$pop, 
    MaxSep = mul.discr.MaxSep(lstX, U)$pop, 
    AverageD = mul.discr.averageD(lstX, U)$pop, 
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop, 
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)
  
})
fiveGumbel_rate <- t(apply(re2, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(fiveGumbel_rate) <- c("Frank", "Gumbel", "Clayton", "Normal")
print(fiveGumbel_rate)
# Clayton
set.seed(seed)
re3 <- replicate(m, expr = {
  
  U <- copula::rCopula(nx, copula = five_Clayton)
  
  frank.param <- coef(copula::fitCopula(copula = frankCopula(dim = p), data = U, method = "mpl"))
  frank.Y <- copula::rCopula(nx, copula = frankCopula(param = frank.param, dim = p))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  
  gumbel.param <-coef(copula::fitCopula(copula = gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.Y <- copula::rCopula(nx, copula = gumbelCopula(param = gumbel.param, dim = p))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  
  
  clayton.param <- coef(copula::fitCopula(copula = claytonCopula(dim = p), data = U, method = "mpl"))
  clayton.Y <- copula::rCopula(nx, copula = claytonCopula(param = clayton.param, dim = p))
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  
  normal.param <- coef(copula::fitCopula(copula = normalCopula(dim= p), data = U, method = "mpl"))
  normal.Y <- copula::rCopula(nx, copula = normalCopula(param = normal.param, dim = p))
  normal.cop <- copula::normalCopula(param = normal.param, dim = p)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, normal.Y)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, normal.cop)
  
  c(discr_linear = mul.discr.linear(lstX, U)$pop, 
    MaxSep = mul.discr.MaxSep(lstX, U)$pop, 
    AverageD = mul.discr.averageD(lstX, U)$pop, 
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop, 
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)
  
})
fiveClayton_rate <- t(apply(re3, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(fiveClayton_rate) <- c("Frank", "Gumbel", "Clayton", "Normal")
print(fiveClayton_rate)
# Normal
set.seed(seed)
re4 <- replicate(m, expr = {
  
  U <- copula::rCopula(nx, copula = five_Normal)
  
  frank.param <- coef(copula::fitCopula(copula = frankCopula(dim = p), data = U, method = "mpl"))
  frank.Y <- copula::rCopula(nx, copula = frankCopula(param = frank.param, dim = p))
  frank.cop <- copula::frankCopula(param = frank.param, dim = p)
  
  gumbel.param <-coef(copula::fitCopula(copula = gumbelCopula(dim = p), data = U, method = "mpl"))
  gumbel.Y <- copula::rCopula(nx, copula = gumbelCopula(param = gumbel.param, dim = p))
  gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
  
  
  clayton.param <- coef(copula::fitCopula(copula = claytonCopula(dim = p), data = U, method = "mpl"))
  clayton.Y <- copula::rCopula(nx, copula = claytonCopula(param = clayton.param, dim = p))
  clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
  
  normal.param <- coef(copula::fitCopula(copula = normalCopula(dim= p), data = U, method = "mpl"))
  normal.Y <- copula::rCopula(nx, copula = normalCopula(param = normal.param, dim = p))
  normal.cop <- copula::normalCopula(param = normal.param, dim = p)
  
  lstX <- list(frank.Y, gumbel.Y, clayton.Y, normal.Y)
  lstCop <- list(frank.cop, gumbel.cop, clayton.cop, normal.cop)
  
  c(discr_linear = mul.discr.linear(lstX, U)$pop, 
    MaxSep = mul.discr.MaxSep(lstX, U)$pop, 
    AverageD = mul.discr.averageD(lstX, U)$pop, 
    Leverage = mul.discr.leverage(lstX, U)$pop,
    Eigen = mul.discr.eigen(lstX, U)$pop, 
    Eigen_reg = mul.discr.eigen.reg(lstX, U, w = rep(1, 3))$pop, 
    Eigen_reg_w = mul.discr.eigen.reg(lstX, U, w = myweight)$pop, 
    # Emp_copulas = mul.Emp.copulas(lstCop, U)$pop,
    NP = mul.NP(lstCop, U)$pop)
  
})

fiveNormal_rate <- t(apply(re4, MARGIN = 1, FUN =  function(arr) {c(mean(arr== "Population 1"), mean(arr=="Population 2"), mean(arr == "Population 3"), mean(arr== "Population 4"))}))
colnames(fiveNormal_rate) <- c("Frank", "Gumbel", "Clayton", "Normal")
print(fiveNormal_rate)
##########################################
##########################################
