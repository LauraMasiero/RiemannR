\name{VectorizeSCM}
\alias{VectorizeSCM}

\title{
Vectorization of SCMs
}
\description{
This function vectorizes sample covariance matrices with the \code{vect} function.
}
\usage{
VectorizeSCM(scm_list)
}

\arguments{
  \item{scm_list}{a list where every element is a list of numeric covariance matrices.
}
}

\value{
It returns a list of vectors, which are the vectorized form of the input sample covariance  matrices.
}

\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\seealso{
\code{\link{vect}},
\code{\link{SCM}}
}

\keyword{ ~SCM }
