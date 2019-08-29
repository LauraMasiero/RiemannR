\name{geometric_mean}
\alias{geometric_mean}

\title{
Computation of SCMs' geometric mean
}
\description{
This function allows to compute the geometric mean of the given sample covariance matrices.
}
\usage{
geometric_mean(scm_list, eps=1/3)
}

\arguments{
  \item{scm_list}{
a list of numeric sample covariance matrices.
}
\item{eps}{parameter in [0,1] that controls the convergence speed. Higher eps values lead to a faster convergence, but the more close to 1 it is, the more irrealistic the underlining hypothesis are (eps=1 indicates data without ANY noise).}
}

\value{
It returns a list containing the geometric mean.
}
\references{
The algorithm used to compute the geometric mean can be found here:
Moakher, Maher (2005). "A differential geometric approach to the geometric mean of symmetric positive-definite matrices". In: SIAM Journal on Ma- trix Analysis and Applications 26.3, pp. 735–747. issn: 08954798. doi:
10.1137/S0895479803436937.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\keyword{ ~mean }
\keyword{ ~geometric}
