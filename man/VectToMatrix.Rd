\name{vectToMatrix}
\alias{vectToMatrix}

\title{
Creation of a vectorized SCMs matrix
}
\description{
This function creates a matrix containing the vectorized form (with the operator \code{vect}) of a list of sample covariance matrices already vectorized.
}
\usage{
vectToMatrix(scm_vect_list)
}

\arguments{
  \item{scm_vect_list}{a list of vectors, which are sample covariance matrices vectorized with the \code{vect} operator.

}
}

\value{
It returns a matrix where the rows are the sample covariance matrices vectorized.
}

\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\seealso{
\code{\link{vect}}
}

\keyword{ ~vect }
\keyword{ ~vectmatr }
