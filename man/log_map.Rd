\name{log_map}
\alias{log_map}

\title{
Logarithmic mapping of SCMs
}
\description{
This function allows to project sample covariance matrices which lies in a Riemannian manifold at the tangent space computed at a given point of the manifold.
}
\usage{
log_map(scm_list, C)
}

\arguments{
  \item{scm_list}{
a list of sample covariance matrices.
}
\item{C}{centroid. The point of the Riemannian manifold in which the tangent space is computed.}
}

\value{
It returns a list containing the matrices projected in the tangent space.
}
\references{
Barachant, Alexandre, Stéphane Bonnet, Marco Congedo e Christian Jutten (2012). "Multiclass brain-computer interface classification by Riemannian geometry". In: IEEE Transactions on Biomedical Engineering. issn: 00189294. doi: 10.1109/TBME.2011.2172210.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\keyword{ ~map }

