## -- filename: CF_bivariate_2.R
library(plotly)
library(copula)
source("Arc_copulas.R")
######################### Miss rate ###########################################
nx <- 500
tau <- 0.7; p <- 2
mycops <- Arc_copula(p = p, tau = tau, df = 2)
bi_Clayton <- mycops$clayton.cop
bi_Frank <- mycops$frank.cop
bi_Gumbel <-  mycops$gumbel.cop
bi_joe <- mycops$joe.cop
bi_t <- mycops$t.cop
bi_Normal <- mycops$normal.cop

set.seed(12)
test1 <- copula::rCopula(nx, copula = bi_Frank)
test2 <- copula::rCopula(nx, copula = bi_Gumbel)
test3 <- copula::rCopula(nx, copula = bi_Clayton)
test4 <- copula::rCopula(nx, copula = bi_joe)
test5 <- copula::rCopula(nx, copula = bi_Normal)
test6 <- copula::rCopula(nx, copula = bi_t)

par(mfrow = c(3, 2), cex.axis = 1.2, 
    cex.lab = 1.2, mar = c(4, 4, 2,1), cex.main = 1.2
)
plot(test1[, 1], test1[, 2], main = "Frank copula", xlab = "", ylab = "")
plot(test2[, 1], test2[, 2], main = "Gumbel copula", xlab = "", ylab = "")
plot(test3[, 1], test3[, 2], main = "Clayton copula", xlab = "", ylab = "")
plot(test4[, 1], test4[, 2], main = "Joe copula", xlab = "", ylab = "")
plot(test5[, 1], test5[, 2], main = "Normal copula", xlab = "", ylab = "")
plot(test6[, 1], test6[, 2], main = "t copula, df = 2", xlab = "", ylab = "")

myBiCopulas <- list("Frank" = bi_Frank, "Gumbel" = bi_Gumbel, 
                    "Clayton" = bi_Clayton, "Joe" = bi_joe, 
                    "Normal" = bi_Normal, "t" = bi_t)
# Plots copulas with different margins 
########################################
library(plotly)
mycolor <- matrix(c(0, "#FFFFFF", 
                    1, "#2171B5"), nrow = 2, byrow = T)
for (k in 1:length(myBiCopulas)){
  cops <- myBiCopulas[[k]]
  print(k)
  dcops <- copula::mvdc(
    cops, margins = c("norm", "norm"),  
    paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
  )
  x1 <- x2 <- seq(-3, 3, length = 100)
  v<-c()
  for (i in x1){
    for (j in x2){
      v<-c(v,copula::dMvdc(c(i,j), dcops) )
    }
  }
  
  f<-t(matrix(v,nrow=100,byrow=TRUE))
  
  g <- plotly::plot_ly(x=x1,y=x2,z=f,type = "contour",colorscale = mycolor, 
                       reversescale = F, showscale = F) %>%
    plotly::layout(xaxis=list(title="x1"),
                   yaxis=list(title="x2"),
                   title=paste(names(myBiCopulas)[k], "copulas"))
  print(g)
}

#######################################
#### Get misrate of Copulas classification using eigenstructure
myBiRate <- list()
length(myBiRate) <- length(myBiCopulas)
names(myBiRate) <- names(myBiCopulas)
######################### Miss rate ###########################################
m <- 2e3; nx <- 500; myweight <- c(2/5, 1/5, 2/5) #change the number
seed <- 1234; 
set.seed(seed)
for (i in 1:length(myBiCopulas)){
  cops <- myBiCopulas[[i]]
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
        copula = copula::normalCopula(dim = p), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::normalCopula(dim = p), data = U, method = "itau"
        )
      })
    )
    normal.cop <- copula::normalCopula(param = normal.param, dim = p)
    normal.Y <- copula::rCopula(nx, copula = normal.cop)
    normal.Y2 <- copula::rCopula(2*nx, copula = normal.cop)
    normal.Y4 <- copula::rCopula(4*nx, copula = normal.cop)
    normal.Y8 <- copula::rCopula(8*nx, copula = normal.cop)
    normal.Y16 <- copula::rCopula(16*nx, copula = normal.cop)
    normal.Y32 <- copula::rCopula(32*nx, copula = normal.cop)
    
    t.param <- coef(tryCatch(
      copula::fitCopula(
        copula = copula::tCopula(dim = p, df.fixed = FALSE), data = U, method = "mpl"
      ), 
      error = function(ex){
        warning("MLE has error - ", ex)
        copula::fitCopula(
          copula = copula::tCopula(dim = p, dispstr = "un"),
          data = U, method = "itau.mpl"
        )
      })
    )
    
    t.cop <- copula::tCopula(
      param = t.param[1], dim = p, df = round(t.param[2])
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
  
  bi_rate_tmp <- t(
    apply(
      re, MARGIN = 1, 
      FUN = function(arr) {
        c(mean(arr== "Population 1"), 
          mean(arr=="Population 2"),
          mean(arr == "Population 3"),
          mean(arr== "Population 4"), 
          mean(arr== "Population 5"), 
          mean(arr== "Population 6")
        )}
    )
  )
  colnames(bi_rate_tmp) <- 
    c("Frank", "Gumbel", "Clayton", "Joe", "Normal", "t")
  myBiRate[[i]] <- bi_rate_tmp
}
tmp <- myBiRate[[1]][, 1]
for (k in 2:6){
  tmp <- cbind(tmp, myBiRate[[k]][, k])
}
tmp
xtable::xtable(tmp, digits = 4)
tmp <- myBiRate[[1]]
# print the results. 
for (k in 1:9){
  # k <- 9
  tmp <- sapply(1:6, function(x) myBiRate[[x]][k, ])
  diag(tmp) <- NA
  tmp1 <- rbind(apply(tmp, MARGIN = 2, function(x) which.max(x)), 
                apply(tmp, MARGIN = 2, function(x) max(x, na.rm = T))
  )
  colnames(tmp1) <- c("Frank", "Gumbel", "Clayton", "Joe", "Normal", "t")
  tmp1[1, ] <- case_when(tmp1[1, ] == 1 ~ "Frank", 
                         tmp1[1, ] == 2 ~ "Gumbel", 
                         tmp1[1, ] == 3 ~ "Clayton",
                         tmp1[1, ] == 4 ~ "Joe", 
                         tmp1[1, ] == 5 ~ "Normal", 
                         tmp1[1, ] == 6 ~ "t", 
  )
  tmp1 <- xtable::xtable(tmp1, digits = 4)
  print(tmp1)
}


