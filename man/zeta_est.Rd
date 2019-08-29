\name{zeta_est}
\alias{zeta_est}

\title{
Estimate of the normalization constant of the Riemann-Gauss distribution
}
\description{
This function computes the normalization constant of the Riemann-Gauss distribution with simulations.
}
\usage{
zeta_est(eta,n_channels,R)
}

\arguments{
  \item{eta}{a vector of numeric values (or a single real value) used to initialize the simulation. \code{eta=-1/2*sigma^2}}
}
\item{n_channels}{a real number indicating the number of channels used in the EEG recordings.}
\item{R}{Number of simulations.}
}

\value{
It returns a numeric vector (if the input \code{eta} was a vector, a number if \code{eta} was a real number) with the estimate of the normalization constant.
}
\references{
Said, Salem, Lionel Bombrun, Yannick Berthoumieu e Jonathan Manton (2015). "Riemannian Gaussian Distributions on the Space of Symmetric Positive Definite Matrices". In: arXiv: 1507 . 01760. url: http://arxiv.org/abs/1507.01760.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}



\examples{
values <- seq(1,2000)^2 #zeta(eta) simulation for a large sample
eta_values <- -10^6/values
zeta <- sapply(eta_values,zeta_est,n_channels=16,10000)
summary(log(zeta[3:2000]))
psi_smooth <- smooth.spline(eta_values[3:2000],log(zeta)[3:2000],nknots=1800) # function smoothing
logzetaprime <- predict(psi_smooth,der=1) # first derivative of psi
}

\keyword{ ~zeta }
\keyword{ ~normalization }
