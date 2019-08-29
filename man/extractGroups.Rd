\name{extractGroups}
\alias{extractGroups}

\title{
Extracting SCMs of the same group from a list
}
\description{
This function allows to extract from a list containing sample covariance matrices the ones of the same group. It takes into consideration only two possible groups, labelled by 0 and 1.
}
\usage{
extractGroups(list, dataset, events=1)
}

\arguments{
  \item{list}{a list of sample covariance matrices to be divided in groups.}
}
\item{dataset}{the dataset of reference. It has to be constructed in a specific manner. See Details below.}
\tem{events}{numeric value: 0 or 1. It is the label of the group to be extracted.}
}
\details{The \code{dataset} must have a specific construction: it must be a list in which every  element is a list itself with two main elements, \code{$epochs} and \code{$labels}.
\code{$epochs} contains tensors with dimensions \code{[n,c,t]}, where n are the number of epochs, c are the number of channels and t is the recording (epoch) length.
\code{$labels} contains the class of the epochs considered (i.e. target or non target).
}
\value{
It returns a list containing only the SCMs that are part of the group identified by \code{events}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}

\seealso{

}

\keyword{ ~groups}
\keyword{ ~extract }
