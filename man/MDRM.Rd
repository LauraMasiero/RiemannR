\name{MDRM}
\alias{MDRM}

\title{
Classification of covariance matrices with MDRM.
}
\description{
This function allows the binary classification of the given covariance matrices with the Minimum distance to Riemannian Mean algorithm.
}
\usage{
MDRM(scm_list, centroids_list)
}

\arguments{
  \item{scm_list}{a list of numeric sample covariance matrices.
}
\item{centroids_list}{a list containing the (geometric) mean of both classes computed on a training set. The first element of the list must be the target class' mean, followed by the mean of the other class.}
}

\value{
It returns a list containing two elements:
\item{classes}{numeric vector cointaining the class indicator.}
\item{distances}{numeric matrix containing the distances among the two centroids, computed with the Riemannian distance.}
}
\references{
Barachant, Alexandre e Marco Congedo (2014). "A Plug&Play P300 BCI Using Information Geometry". In: url: http://arxiv.org/abs/1409.0107.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}



\seealso{
 \code{\link{RaoDist}}
}

\keyword{ ~Min }
\keyword{ ~Distance }
