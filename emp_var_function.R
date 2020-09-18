emp_var_function <- function(mu, scm_list) {
  P <- length(scm_list)
  P^(-1)*Reduce("+",lapply(scm_list, function(x) RaoDist(x, mu)))
}

eta <- function(sigma) {
  -1/(2*sigma^2)
} 

