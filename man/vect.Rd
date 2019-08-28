\name{vect}
\alias{vect}
\title{
Vectorization of covariance matrices
}
\description{
This function vectorizes sample covariance matrices with the formulation presented by
}
\usage{
vect(x)
}

\arguments{
  \item{x}{
a numeric matrix.
}
}
\details{
This type of vectorization takes the upper triangular part of a matrix, multiplying the non-diagonal elements with an appropriate weight. For more details see References.
}
\value{
It returns the vectorized form of the input matrix.
}
\references{
Barachant, Alexandre, Stéphane Bonnet, Marco Congedo e Christian Jutten (2013). "Classification of covariance matrices using a Riemannian-based kernel for BCI applications". In: Neurocomputing. issn: 09252312. doi: 10.1016/j.neucom.2012.12.039.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\keyword{ ~vec }

