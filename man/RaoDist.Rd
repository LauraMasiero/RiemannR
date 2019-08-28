\name{RaoDist}
\alias{RaoDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Riemannian Distance
}
\description{
This function computes the Riemannian distance between two matrices with the Rao-Fisher metric.
}
\usage{
RaoDist(m1,m2)
}

\arguments{
  \item{m1, m2}{numeric matrices.
}
}
\details{
The distance is computed with the use of the spectral deccomposition.
}
\value{
It returns a numeric value, the distance.
}
\references{
Moakher, Maher (2005). "A differential geometric approach to the geometric mean of symmetric positive-definite matrices". In: SIAM Journal on Matrix Analysis and Applications 26.3, pp. 735–747. doi: 10.1137/S0895479803436937.

Said, Salem, Lionel Bombrun, Yannick Berthoumieu e Jonathan Manton (2015). "Riemannian Gaussian Distributions on the Space of Symmetric Positive Definite Matrices". In: arXiv: 1507 . 01760. url: http://arxiv.org/abs/1507.01760.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\keyword{ ~Rao }
\keyword{ ~Dist }
