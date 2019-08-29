\name{RG_distr}
\alias{RG_distr}

\title{
Riemann-Gauss density function
}
\description{
This function computes the Riemann-Gauss density on a given value.
}
\usage{
RG_distr(mu,sigma, p, scm_matrix)
}

\arguments{
  \item{mu}{a numeric matrix indicating the distribution's mean.
}
\item{sigma}{a numeric value indicating the dispersion of the distribution.}
\item{p}{a numeric value representing the number of channels used in the EEG recording.}
\item{scm_matrix}{a matrix in which the density value will be computed.}
}

\value{
It returns a numeric value indicating the estimate density in \code{scm_matrix}.
}
\references{
Said, Salem, Lionel Bombrun, Yannick Berthoumieu e Jonathan Manton (2015). "Riemannian Gaussian Distributions on the Space of Symmetric Positive Definite Matrices". In: arXiv: 1507 . 01760. url: http://arxiv.org/abs/1507.01760.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\keyword{ ~Riemann }
\keyword{ ~Gauss }
\keyword{ ~RG }
\keyword{ ~distribution }
