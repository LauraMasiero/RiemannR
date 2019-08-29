\name{scm_transform}
\alias{scm_transform}

\title{SCMs transformation to use a linear SVM

}
\description{
This function allows to transform the given sample covariance matrices in a matrix usable in a SVM with linear kernel. This is a different way to implement a Riemannian kernel SVM for SCMs classification.
}
\usage{
scm_transform(scm_list, centroid)
}

\arguments{
  \item{scm_list}{a list of numeric sample covariance matrices.
}
\item{centroid}{a numeric matrix used as reference for the tangent space mapping. See Details below.}
}
\details{
The SCMs transformation is achieved with a logaritmic mapping of the matrices, that's to say the matrices are projected from the Riemannian manifold in the tangent euclidean space. The tangent space is computed in a reference point: the \code{centroid}. For the formulation details see the  References below.
}
\value{
It returns a numeric matrix conatining all the upper triangular parts of the SCM in input, vectorized with the \code{vect} operator.
}
\references{
Barachant, Alexandre, Stéphane Bonnet, Marco Congedo e Christian Jutten, (2013). "Classification of covariance matrices using a Riemannian-based kernel for BCI applications". In: Neurocomputing. issn: 09252312. doi: 10.1016/j.neucom.2012.12.039.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\seealso{
\code{\link{vect}}
}

\keyword{ ~transform }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~vect }% __ONLY ONE__ keyword per line
