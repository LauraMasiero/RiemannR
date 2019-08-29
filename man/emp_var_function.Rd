\name{emp_var_function}
\alias{emp_var_function}

\title{
Empirical variance function
}
\description{
This function computes the empirical variance function as described in Said et al., 2015.
}
\usage{
emp_var_function(mu, scm_list)
}

\arguments{
  \item{mu}{the geometric mean of the sample covariance matrices in \code{scm_list}. It is in fact the mean of the Riemann-Gauss distribution described in Said et al., 2015. For more details see the References below.
}
\item{scm_list}{a list of numeric sample covariance matrices.}
}


\value{
It returns a numeric value, that is the empirical variance function estimated.
}
\references{
Said, Salem, Lionel Bombrun, Yannick Berthoumieu e Jonathan Manton (2015). "Riemannian Gaussian Distributions on the Space of Symmetric Positive Definite Matrices". In: arXiv: 1507 . 01760. url: http://arxiv.org/abs/1507.01760.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\seealso{
\code{\link{RG_distr}},
}

\keyword{ ~emp }
\keyword{ ~var }
