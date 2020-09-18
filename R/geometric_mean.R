geometric_mean <- function(scm_list, eps=1/3){
  P <- length(scm_list)
  cat("Initializing...", "\n")
  geo_mean <- P^(-1)*Reduce("+", scm_list) #inizializzazione
  cat("Estimating geometric mean...", "\n")
  repeat{
    x <- geo_mean
    S_hat <-P^(-1)*Reduce("+", log_map(scm_list, x)) #media aritmetica nello spazio tangente
    geo_mean <- exp_map(S_hat, x)
    if(norm(S_hat, type="F")<eps) break
  }
  cat("Geometric mean estimated.", "\n")
  return(geo_mean)
}