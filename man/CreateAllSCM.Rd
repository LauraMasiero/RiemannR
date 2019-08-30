\name{createAllSCM}
\alias{createAllSCM}

\title{
Creation of sample covariance matrices
}
\description{
This function creates the sample covariance matrices of multiple matrices created with or without the function \code{SuperTrial}.
}
\usage{
createAllSCM(list, n=8, st=T, tidy=FALSE, col=17)
}

\arguments{
  \item{list}{a list of matrices with a specific construction (if \code{st=F}) or a list of matrices constructed by the function \code{SuperTrial} (if \code{st=T}). See Details below.
}
\item{n}{a numeric value indicating how many SCM have to be computed.}
\item{st}{boolean. If \code{st=T} the elements of the \code{list} in input are matrices constructed with the function \code{SuperTrial}, if \code{st=F} the matrices have a particular construction, described in Details below. }
\item{tidy}{option to remove specific channels. If \code{tidy=TRUE} no channels will be removed to construct the Super trial, if \code{tidy=FALSE} then the columns of the channels to be removed have to be specified with the option \code{col}. }
\item{col}{columns to be removed (when \code{tidy=FALSE}).}
}


\details{The initial list must have a specific construction: every list element has to be a list itself with two main elements, \code{$epochs} and \code{$labels}.
\code{$epochs} contains tensors with dimensions \code{[n,c,t]}, where n are the number of epochs, c are the number of channels and t is the recording (epoch) length.
\code{$labels} contains the class of the epochs considered (i.e. target or non target).

}
\value{
It returns a list containing the sample covariance matrices computed.
}
\references{
Barachant, Alexandre, Stéphane Bonnet, Marco Congedo e Christian Jutten (2012). "Multiclass brain-computer interface classification by Riemannian geometry". In: IEEE Transactions on Biomedical Engineering. issn: 00189294. doi: 10.1109/TBME.2011.2172210.

Barachant, Alexandre e Marco Congedo (2014). “A Plug&Play P300 BCI Using Information Geometry”. url: http://arxiv.org/abs/1409.0107.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\seealso{
\code{\link{SuperTrial}}
}

\keyword{ ~SCM }

