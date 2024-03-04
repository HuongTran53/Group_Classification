## -- Filename: discr_NP.R
# Perform optimal classification rule, Neyman-Pearson
NP <- function(U, lst_f, use.log = T){
  re <- lapply(lst_f, function(f) {
    ifelse(
      use.log == T, 
      sum(log(apply(U, MARGIN = 1, f))), 
      sum(apply(U, MARGIN = 1, f))
    )
  })
  return(which.max(re))
}