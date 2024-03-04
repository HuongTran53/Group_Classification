## --Filename: Arc_copulas.R
## Generate some copulas with given tau correlation and dimension 
Arc_copula <- function(p, tau, df = 4){
  require(copula)
  normal.cor <- sin(pi*tau/2)
  clayton.alpha <- (2 *tau)/(1- tau)
  gumbel.alpha <- 1 /(1 - tau)
  ####
  fun <- function(theta){
    1/theta * integrate(
      function(t) t/(exp(t) - 1), lower = 0, upper = theta
      )$value
  }
  frank.alpha  <- uniroot(function(theta) 1 + 4/theta *(fun(theta) - 1) - tau, 
                          lower = -100, upper = 100)$root
  
  fun.joe <- function(theta){
    1 + 4/theta^2 * integrate(
      function(t) t*log(t)*(1 - t)^{2*(1 - theta)/theta}, lower = 0, upper = 1
      )$value
  }
  joe.alpha <- uniroot(
    function(theta) tau - fun.joe(theta), lower = 1, upper = 100
    )$root
  
  frank.cop <- copula::frankCopula(param = frank.alpha, dim = p)
  gumbel.cop <- copula::gumbelCopula(gumbel.alpha, dim = p)
  clayton.cop <- copula::claytonCopula(param = clayton.alpha ,dim = p)
  joe.cop <- copula::joeCopula(param = joe.alpha, dim = p)
  normal.cop <- copula::normalCopula(param = normal.cor, dim = p)
  t.cop <- copula::tCopula(param = normal.cor, dim = p, df = df)
  
  return(list(frank.cop = frank.cop, 
              gumbel.cop = gumbel.cop, 
              clayton.cop = clayton.cop, 
              joe.cop = joe.cop,
              normal.cop = normal.cop, 
              t.cop = t.cop))
}

