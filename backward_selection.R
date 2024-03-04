## --Filename: backward_selection.R
backward_selection <- function(X1, X2){
  p <- ncol(X1)
  nX1 <- nrow(X1)
  nX2 <- nrow(X2)
  
  capH1 <- (nX1/nX2) *solve(t(X1) %*% X1) %*% (t(X2) %*% X2)
  detcapH1 <- det(capH1); trcapH1 <- sum(diag(capH1))
  zeta1 <- detcapH1^(1/p)/(trcapH1/p)
  capH2 <- (nX2/nX1) *solve(t(X2) %*% X2) %*% (t(X1) %*% X1)

  detcapH2 <- 1/detcapH1; 
  trcapH2 <- sum(diag(capH2))
  zeta2 <- detcapH2^(1/p)/(trcapH2/p)
  
  zeta <- sqrt(prod(zeta1, zeta2))
  zetamin <- 0 
  allcol <- 1:p
  capY1 <- X1
  capY2 <- X2
  q <- p
  rm <- c() # keep track of removed columns
  
  while (zetamin < zeta){
    lstzeta <- c()
    for (j in 1:q){
      Y1 <- capY1[, -j]
      Y2 <- capY2[, -j]
      y1 <- capY1[, j]
      y2 <- capY2[, j]
      
      H1 <- (nX1/nX2) * solve(t(Y1) %*% Y1) %*% (t(Y2) %*% Y2)
      H2 <- (nX2/nX1) * solve(t(Y2) %*% Y2) %*% (t(Y1) %*% Y1)
      detH1 <- det(H1)
      detH2 <- det(H2)
      trH1 <- sum(diag(H1))
      trH2 <- sum(diag(H2))
      
      zetaH1 <- (q - 1) * detH1^(1/(q - 1))/trH1
      zetaH2 <- (q - 1) * detH2^(1/(q - 1))/trH2
    
      myzeta <- sqrt(prod(zetaH1, zetaH2))
      lstzeta <- c(lstzeta, myzeta)
    }
    zetamin <- min(lstzeta)
    if (zetamin < zeta){
      index <- which.min(lstzeta)             
      rm <- c(rm, allcol[index])
      allcol <- allcol[-index]
      q <- q - 1
      zeta <- zetamin
      zetamin <- 0
      capY1 <- capY1[, -index]
      capY2 <- capY2[, -index]
    } else {
      break
    }
  }
  return(list(namecol = allcol, zeta = zeta, remove.col = rm))
}
