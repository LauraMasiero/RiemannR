\name{LR_classification}
\alias{LR_classification}

\title{
Binary classification with RG ditribution
}
\description{
This function computes binary classification using a likelihood ratio, based on Riemann-Gauss distribution.
}
\usage{
LR_classification(mu, sigma, scm_list)
}

\arguments{
  \item{mu}{a list containing the (geometric) mean of the two different classes.
}
\item{sigma}{a numeric vector containing the dispersion parameter of the two classes.}
\item{scm_list}{a list containoing the sample covariance matrices to classify.}
}

\value{
It returns a list containing two elements:
\item{lik_ratios}{a list where each element is the likelihood ratio}
\item{classes}{a numeric vector indicating the class estimated for every scm matrix.}
}


\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\seealso{
\code{\link{RG_distr}
}

\keyword{ ~LR }
\keyword{ ~classification }
