##--Filename: CF_varsel.R
### Sampe covariance, diff mean 
mu <- c(rep(-10, 3), rep(-5, 2), rep(0, 10), rep(7, 3), rep(15, 2))
p <- length(mu)
set.seed(123)
mu1 <- sample(mu)
mu2 <- sample(mu)
plot(abs(mu1 - mu2))
text(1:p, abs(mu1 - mu2) , as.character(1:p), pos=1)

set.seed(1111)
X1 <- MASS::mvrnorm(500, mu1, diag(p))
X2 <- MASS::mvrnorm(700, mu2, diag(p))

forward_selection(X1, X2)
backward_selection(X1, X2)

### Same mean, diff covariance matrix
n1 <- 500; n2 <- 700
q1 <- 10
S <- diag(q1);
eigen(S, only.values = T)
set.seed(1234)
Y1 <- MASS::mvrnorm(n1, rep(0, q1), S)
Y2 <- MASS::mvrnorm(n2, rep(0, q1), S)
###
q2 <- p - q1
###
set.seed(123)
J <- c(0, 1, 2)
A <- matrix(sample(J, size = q2*q2, replace = T), nrow = q2, byrow = T)
diag(A)[diag(A) == 0] <- rep(1, sum(diag(A)==0))
S1 <- t(A) %*% A
set.seed(123)
X1 <- MASS::mvrnorm(n1, rep(0, q2), S1)
####
J <- c(-2, -1, 0, 1, 2)
set.seed(123)
A <- matrix(sample(J, size = q2*q2, replace = T), nrow = q1, byrow = T)
diag(A)[diag(A) <= 0] <- rep(1, sum(diag(A)==0))
S2 <- t(A) %*% A
set.seed(123)
X2 <- MASS::mvrnorm(n2, rep(0, q2), S2)

X1 <- cbind(Y1, X1)
X2 <- cbind(Y2, X2)
forward_selection(X1, X2)
backward_selection(X1, X2)

