\name{exp_map}
\alias{exp_map}

\title{
Exponential mapping of SCMs
}
\description{
This function allows to project bacj to the Riemannian manifold sample covariance matrices which were mapped in tangent space computed at a given point of the manifold.
}
\usage{
exp_map(mapped_values, C)
}
\arguments{
  \item{mapped_values}{
a list of sample covariance matrices previously mapped in the tangent space computed at C.
}
\item{C}{centroid. The point of the Riemannian manifold in which the tangent space was computed.}
}
\value{
It returns a list containing the matrices projected in the Riemannian manifold.
}
\references{
Barachant, Alexandre, Stéphane Bonnet, Marco Congedo e Christian Jutten (2012). "Multiclass brain-computer interface classification by Riemannian geometry". In: IEEE Transactions on Biomedical Engineering. issn: 00189294. doi: 10.1109/TBME.2011.2172210.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\keyword{ ~map }
