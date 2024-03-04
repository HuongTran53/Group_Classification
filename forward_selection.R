## --filename: forward_selection.R
forward_selection <- function(X1, X2, mincol = 2){
  # mincol: Number of dersied columns
  # choose the best columns at every step to add on
  p <- ncol(X1)
  nX1 <- nrow(X1)
  nX2 <- nrow(X2)
  q <- 2
  col <- combn(1:p, q)
  
  sX1 <- 1/sqrt(nX1) * X1
  sX2 <- 1/sqrt(nX2) * X2
  
  # Initilization
  mycol <- col[, 1]
  H1 <- solve(
    t(sX1[, mycol]) %*% sX1[, mycol]
  ) %*% (
    t(sX2[, mycol]) %*% sX2[, mycol]
  )
  detH1 <- det(H1); trH1 <- sum(diag(H1))
  zeta1 <- detH1^(1/q)/(trH1/q)
  H2 <- solve(
    t(sX2[, mycol]) %*% sX2[, mycol]
  ) %*% (
    t(sX1[, mycol]) %*% sX1[, mycol]
  )
  # detH2 <- det(H2);
  detH2 <- 1/detH1; 
  trH2 <- sum(diag(H2))
  zeta2 <- detH2^(1/q)/(trH2/q)
  # zeta <- min(zeta1, zeta2)
  zeta <- sqrt(prod(zeta1, zeta2))
  
  for (subcol in 2:ncol(col)){
    mycol_tmp <- col[, subcol]
    H1_tmp <- solve(
      t(sX1[, mycol_tmp]) %*% sX1[, mycol_tmp]
    ) %*% (
      t(sX2[, mycol_tmp]) %*% sX2[, mycol_tmp]
    )
    detH1_tmp <- det(H1_tmp); trH1_tmp <- sum(diag(H1_tmp))
    zeta1_tmp <- detH1_tmp^(1/q)/(trH1_tmp/q)
    
    H2_tmp <- solve(
      t(sX2[, mycol_tmp]) %*% sX2[, mycol_tmp]
    ) %*% (
      t(sX1[, mycol_tmp]) %*% sX1[, mycol_tmp]
    )
    detH2_tmp <- 1/detH1_tmp
    trH2_tmp <- sum(diag(H2_tmp))
    zeta2_tmp <- detH2_tmp^(1/q)/(trH2_tmp/q)
    zeta_tmp <- sqrt(prod(zeta1_tmp, zeta2_tmp))
    
    if (zeta_tmp < zeta){
      mycol <- mycol_tmp
      H1 <- H1_tmp
      H2 <- H2_tmp
      detH1 <- detH1_tmp; trH1 <- trH1_tmp
      detH2 <- detH2_tmp; trH2 <- trH2_tmp
      zeta <- zeta_tmp
    }
  }
  
  Y1 <- sX1[, mycol]
  Y2 <- sX2[, mycol]
  othercol <- c(1:p)[-mycol]
  # q <- q; 
  zetamin <- 0
  # Searching new columns to append 
  while (zetamin < zeta | q < mincol){
    #count <- count + 1 # Avoid the case q is best, 
    # then algorithm will stop after searching in p - q
    Z1 <- Y1 %*% solve(t(Y1) %*% Y1) %*% t(Y1)
    Z2 <- Y2 %*% solve(t(Y2) %*% Y2) %*% t(Y2)
    iZ1 <- diag(nX1) - Z1
    iZ2 <- diag(nX2) - Z2
    Z12 <- Y1 %*% solve(t(Y1) %*% Y1) %*% t(Y2)
    Z21 <- Y2 %*% solve(t(Y2) %*% Y2) %*% t(Y1)
    lstzeta <- c()
    lstdet1 <- lstdet2 <- c()
    lsttr1 <- lsttr2 <- c()
    for (j in othercol){
      ####### (sX1sX2)-1(sX2sX2)
      B <- t(sX1[,j]) %*% iZ1 %*% sX1[,j]
      ZyyZ <- t(Z12) %*% sX1[,j] %*% t(sX1[,j]) %*% Z12
      trcapH1 <- trH1 + 1/B * (
        sum(diag((ZyyZ))) - 
          2 * t(sX2[,j]) %*% t(Z12) %*% sX1[,j] +
          t(sX2[,j]) %*% sX2[,j]
      )
      detcapH1 <- detH1 * (
        t(sX2[,j]) %*% iZ2 %*% sX2[, j]
      )/(
        t(sX1[,j]) %*% iZ1 %*% sX1[, j]
      )
      zetacapH1 <- (q + 1) * detcapH1^(1/(q + 1))/trcapH1
      lstdet1 <- c(lstdet1, detcapH1)
      lsttr1 <- c(lsttr1, trcapH1)
      ####### (sX2X2)-1(x1X1)
      B <- t(sX2[,j]) %*% iZ2 %*% sX2[,j]
      ZyyZ <- t(Z21) %*% sX2[,j] %*% t(sX2[,j]) %*% Z21
      trcapH2 <- trH2 + 1/B * (
        sum(diag((ZyyZ))) - 
          2 * t(sX1[,j]) %*% t(Z21) %*% sX2[,j] + t(sX1[,j]) %*% sX1[,j]
      )
      detcapH2 <- 1/detcapH1
      zetacapH2 <- (q + 1) * detcapH2^(1/(q + 1))/trcapH2
      lstdet2 <- c(lstdet2, detcapH2)
      lsttr2 <- c(lsttr2, trcapH2)
      
      myzeta <- sqrt(prod(zetacapH1, zetacapH2))
      lstzeta <- c(lstzeta, myzeta)
    }
    zetamin <- min(lstzeta)
    if (zetamin < zeta | q < mincol){
      index <- which.min(lstzeta)
      mycol <- c(mycol, othercol[index])
      othercol <- c(1:p)[-mycol]
      zeta <- zetamin
      zetamin <- 0
      Y1 <- sX1[, mycol]
      Y2 <- sX2[, mycol]
      detH1 <- lstdet1[index]
      trH1 <- lsttr1[index]
      
      detH2 <- lstdet2[index]
      trH2 <- lsttr2[index]
      q <- q + 1; 
    } else {
      break
    }
  }
  return(list(namecol = mycol, zeta = zeta))
}
