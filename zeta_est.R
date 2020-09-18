eta <- function(sigma) -1/(2*sigma^2)

zeta_est <- function(eta,n_channels,R){
  require(mvtnorm)
  mu <- rep(0,n_channels)
  Sigma <- diag(-1/eta,n_channels)
  set.seed(123)
  val <- rmvnorm(R,mu,Sigma)
  temp <- rep(1,R)
  for (i in 1:(n_channels-1))
    for (j in (i+1):n_channels)
      temp <- temp*(sinh(abs(val[,i]-val[,j])/2))
  omegap <- ((((pi^((n_channels^2)/2))*8^((n_channels*(n_channels-1))/4)))/(factorial(n_channels)*lmvgamma(n_channels/2,n_channels)))
  res <- omegap*mean(temp)
  res <- res*(2*pi)^(n_channels/2)*(det(Sigma))^(1/2)
  return(res)
}