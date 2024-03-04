## -- Filename: Data_water.R
water <- 
  read.csv("/Users/huongtran/OU /Dissertation/Datasets/water_potability.csv")
library(dplyr)
water <- na.omit(water)
X1 <- water %>% dplyr::filter(Potability == 1) %>% select(-c("Potability"))
X0 <- water %>% dplyr::filter(Potability == 0) %>% select(-c("Potability"))
S1 <- cov(X1)
S0 <- cov(X0)
# S_compare <- (S1 - S0)/S1
S_compare <- S0/S1
par(cex.axis = 1.2, cex.lab = 1.2, mar = c(4, 4.1, 2,1), cex.main = 1.2)
plot(c(S_compare[upper.tri(S_compare, diag = T)]), 
     ylab = "Ratio between elements in sample covariance matrices",
     lwd = 2, pch = 16, col = "#2171b5", yaxt="n",
     xlab = "Index of distinct elements in covariance matrix"
     )
axis(2, at = c(-15, -10, -1.5, 0, 1.5, 5, 10, 15), las=2)
abline(h = 1.5, lty = "dotted", col = "orange", lwd = 1.5)
abline(h = -1.5, lty = "dotted", col = "orange", lwd = 1.5)

e1 <- eigen(cov(X1), only.values = T)$values
e0 <- eigen(cov(X0), only.values = T)$values
par(cex.axis = 1.2, cex.lab = 1.2, mar = c(4, 4.1, 2,1), cex.main = 1.2)
plot(e0/e1, ylim = c(0.5, 1.5), 
     ylab = "Ratio between eigenvalues of sample covariance matrices", 
     lwd = 2, pch = 16, col = "#2171b5", yaxt="n", 
     xlab = "Index of eigenvalues")
axis(2, at = c( .5, .7, 1, 1.3, 1.5), las=2)

mu1 <- colMeans(X1); mu0 <- colMeans(X0)
par(cex.axis = 1.2, cex.lab = 1.2, mar = c(4, 4.1, 2,1), cex.main = 1.2)
plot(mu0/mu1, ylim = c(0.5, 1.5), 
     ylab = "Ratio between elements in sample means", 
     lwd = 2, pch = 16, col = "#2171b5", yaxt="n", 
     xlab = "Index of elements in sample means")
axis(2, at = c( .5, .7, 1, 1.3, 1.5), las=2)
##############################
m <- 2e3
# m <- 1e2
set.seed(1234)
re <- replicate(m, expr ={
  split1 <- rbinom(nrow(X1), 1, 0.3)
  sX1 <- as.matrix(X1[split1 == 0, ]) #Labeled sample
  U1 <- as.matrix(X1[split1 == 1,]) # Unlabeled sample 
  
  split0 <- rbinom(nrow(X0), 1, 0.3)
  sX0 <- as.matrix(X0[split0 == 0,])
  U0 <- as.matrix(X0[split0 == 1, ])
  
  re.linear1 <- discr.linear(X1 = sX1, X2 = sX0, U = U1, true.params = F)$pop 
  re.linear2 <- discr.linear(X1 = sX1, X2 = sX0, U = U0, true.params = F)$pop 
  
  re.MaxSep1 <- discr.MaxSep(X1 = sX1, X2 = sX0, U = U1, true.params = F)$pop
  re.MaxSep2 <- discr.MaxSep(X1 = sX1, X2 = sX0, U = U0, true.params = F)$pop
  
  re.avgD1 <- discr.averageD(X1 = sX1, X2 = sX0, U = U1)$pop 
  re.avgD2 <- discr.averageD(X1 = sX1, X2 = sX0, U = U0)$pop 
  
  re.leverage1 <- discr.leverage(X1 = sX1, X2 = sX0, U = U1)$pop
  re.leverage2 <- discr.leverage(X1 = sX1, X2 = sX0, U = U0)$pop
  
  discr_type_1_1 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U1, type =1)$pop 
  discr_type_1_2 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U0, type =1)$pop 
  
  discr_type_2_1 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U1, type =2)$pop 
  discr_type_2_2 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U0, type =2)$pop 
  
  discr_type_3_1 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U1, type =3)$pop 
  discr_type_3_2 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U0, type =3)$pop 
  
  discr_type_4_1 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U1, type =4)$pop 
  discr_type_4_2 <- discr.antieigen(X1 = sX1, X2 = sX0, U = U0, type =4)$pop 
  
  tmp <- c(
    linear11 = re.linear1,
    linear22 = re.linear2,
    MaxSep11 = re.MaxSep1,
    MaxSep22 = re.MaxSep2,
    avgD11 = re.avgD1, 
    avgD22 = re.avgD2, 
    leverage11 = re.leverage1, 
    leverage22 = re.leverage2, 
    ###########################
    type1_11 = discr_type_1_1, 
    type1_22 = discr_type_1_2,
    type2_11 = discr_type_2_1, 
    type2_22 = discr_type_2_2,
    type3_11 = discr_type_3_1, 
    type3_22 = discr_type_3_2,
    type4_11 = discr_type_4_1, 
    type4_22 = discr_type_4_2
  )
  tmp[tmp == "Can not decide"] <- sample(c("Population 1", "Population 2"), 
                                         size = sum(tmp == "Can not decide"), 
                                         replace = T)
  tmp
})
re_Water <- rbind(apply(re == "Population 1", 1, mean), 
             apply(re == "Can not decide", 1, mean), 
             apply(re == "Population 2", 1, mean))
###########################
### Variable Selection ###########
X1 <- as.matrix(X1)
X0 <- as.matrix(X0)
forward_selection(X1, X0)
backward_selection(X1, X0)
