\name{SuperTrial}
\alias{SuperTrial}

\title{
Super Trial
}
\description{
This function creates a super trial (Barachant and Congedo, 2014) to use in ERP signal analysis for BCI.
}
\usage{
SuperTrial(list, n=16, length_epochs=129, tidy=F, col=17)
}

\arguments{
  \item{list}{a list of matrices with a specific construction. See Details below.
}
\item{n}{number of channels used to record the EEG signal.}
\item{length_epochs}{number of time points of the registered signal (its length).}
\item{tidy}{option to remove specific channels. If \code{tidy=TRUE} no channels will be removed to construct the Super trial, if \code{tidy=FALSE} then the columns of the channels to be removed have to be specified with the option \code{col}. }
\item{col}{columns to be removed (when \code{tidy=FALSE}).}
}


\details{The initial list must have a specific construction: every list element has to be a list itself with two main elements, \code{$epochs} and \code{$labels}.
\code{$epochs} contains tensors with dimensions \code{[n,c,t]}, where n are the number of epochs, c are the number of channels and t is the recording (epoch) length.
\code{$labels} contains the class of the epochs considered (i.e. target or non target).} 

\value{
It returns a list:
 \item{mean_epoch}{a matrix that is the signal mean of every channel considered.}
 \item{super_trial}{a list with the super trials constructed from the initial list.}
} 

\references{
Barachant, Alexandre e Marco Congedo (2014). “A Plug&Play P300 BCI Using Information Geometry”. url: http://arxiv.org/abs/1409.0107.
}
\author{
Laura Masiero, email: laura.masiero.10@gmail.com
}


\keyword{ ~Super }
\keyword{ ~Trial }
